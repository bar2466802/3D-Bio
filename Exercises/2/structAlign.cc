/**
#################################################################
# FILE : structAlign.cc
# WRITER : Bar Melinarskiy & Rachel Ben Hamozeg
# EXERCISE : 3D-Bio ex2 2021
# DESCRIPTION: A program that detects an (almost) largest structural alignment between
two proteins or two RNA molecules.
#################################################################
 */

#include "Vector3.h"
#include "Triangle.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include <chrono>
#include <iostream>

//Constants
#define NUMBER_OF_ARGS 4
#define EPSILON_INDEX 1
#define TARGET_INDEX 2
#define MODEL_INDEX 3

// functions
/**
 * Move the given molecule to the center of mass
 * @param molModel the atoms of the given molecule
 */
Vector3 moveToCenter(Molecule<Atom> *molModel)
{
	// calculate center of mass
	Vector3 vectModelMass(0, 0, 0);
	for (unsigned int i = 0; i < molModel->size(); i++)
	{
		vectModelMass += (*molModel)[i].position();
	}
	vectModelMass /= molModel->size();

	// transform the molecules to the center of the coordinate system
	(*molModel) += (-vectModelMass);

	return vectModelMass;
}

/**
 * Create the updated mol model with the given transformation
 * @param sourceFile the PDB file of the protein / RNA
 * @param rtrans the transformation to apply
 * @return the molecule model after applying the transformation
 */
Molecule<Atom> createUpdatedMolModel(std::ifstream *sourceFile, RigidTrans3 rtrans)
{
	//Read the PDB file with all the molecules
	Molecule<Atom> molModelForExport;
	molModelForExport.readPDBfile(*sourceFile, PDB::AllSelector());

	// calculate center of mass
	moveToCenter(&molModelForExport);

	//Preform rotations
	molModelForExport *= rtrans;

	return molModelForExport;
}

/**
* Create the PDB file with the updated model
* @param fileName export file name
* @param molModelToExport mol model to export to the PDB file
*/
void exportToPDB(std::string fileName, Molecule<Atom> molModelToExport)
{
	std::ofstream of(fileName);
	if (of.is_open())
	{
		of << molModelToExport << std::endl;
		of.flush();
		of.close();
		std::cout << fileName + " was created successfully!" << std::endl;
	}
	else
	{
		std::cerr << "Failed to open file " + fileName << std::endl;
		exit(-1);
	}
}

/**
* Check if the given PDB file holds an RNA molecule
* @param fileName given PDB file to test
* @return true if this is indeed an RNA molecule, false otherwise.
*/
bool checkIfRNA(std::string fileName)
{
	const std::string RNA = " RNA ";
	std::ifstream file(fileName); //open the file
	if (file.is_open()) //check the file was opened successfully
	{
		//Read only the first line - the header
		std::string line;
		getline(file, line);
		return line.find(RNA, 0) != std::string::npos;
	}
	else // an error occurred while opening the given file
	{
		std::cerr << "Failed to open file " + fileName << std::endl;
		exit(-1);
	}
}

/**
* Open the given PDB file
* @param fileName given PDB file to open
* @return stream of the given PDB file.
*/
std::ifstream openPDBfile(std::string fileName)
{
	std::ifstream file(fileName);
	if (!file)
	{
		std::cout << "File " << fileName << " does not exist." << std::endl;
		exit(0);
	}
	else
	{
		std::cout << "File " << fileName << " was successfully found." << std::endl;
	}

	return file;
}

int main(int argc, char *argv[])
{
	// measure the run time
	auto start = std::chrono::system_clock::now();

	if (argc != NUMBER_OF_ARGS)
	{
		std::cerr << "Usage: " << argv[0] << "epsilon target.pdb model.pdb " << std::endl;
		exit(1);
	}

	//********Parameters********************
	float epsilon = atof(argv[EPSILON_INDEX]); // distance threshold on atoms in correspondence
	if (epsilon < 0)
	{
		std::cerr << "Usage: " << argv[0] << "epsilon must be a positive float " << std::endl;
		exit(1);
	}

	std::cout << "Given epsilon: " << epsilon << std::endl;

	// read the two files into Molecule
	Molecule<Atom> molModel, molTarget;

	std::ifstream fileTarget = openPDBfile(argv[TARGET_INDEX]);
	std::ifstream fileModel = openPDBfile(argv[MODEL_INDEX]);

	if (checkIfRNA(argv[MODEL_INDEX]))
	{
		//This is indeed an RNA molecule so read only the phosphate atoms
		molTarget.readPDBfile(fileTarget, PDB::PSelector());
		molModel.readPDBfile(fileModel, PDB::PSelector());
	}
	else
	{
		//This is not an RNA molecule so read only the Cα atoms
		molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
		molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
	}

	//Check the models are not empty
	if (molModel.empty() || molTarget.empty())
	{
		std::cout << "No molecules were found in PDB file" << std::endl;
		return -1;
	}

	// calculate center of mass for both molecules and moved them there
	Vector3 vectTargetMass = moveToCenter(&molTarget);
	Vector3 vectModelMass = moveToCenter(&molModel);


	// next we insert the target molecule into hash
	// this will help us to find atoms that are close faster
	GeomHash<Vector3, int> gHash(3,
								 epsilon); // 3 is a dimension and m_fDistThr is the size of the hash
	// cube
	for (unsigned int i = 0; i < molTarget.size(); i++)
	{
		// coordinate is the key to the hash, we store atom index
		gHash.insert(molTarget[i].position(), i);
	}

	// now we define transformations using triangles formed by 3 consecutive Cα atoms
	//for proteins or phosphate atoms for nucleic acids. Over all |pdb1|x|pdb2| transformations.
	// We will choose the best alignment given from there transformations.
	unsigned int iMaxSize = 0;
	Match bestMatch;
	RigidTrans3 rtransBest;
	//
	for (unsigned int indexTarget = 0; indexTarget < molTarget.size() - 3; indexTarget++)
	{
		//Generate triangles
		Triangle *triangleTarget = new Triangle(molTarget[indexTarget],
												molTarget[indexTarget + 1],
												molTarget[indexTarget + 2]);

		for (unsigned int indexModel = 0; indexModel < molModel.size() - 3; indexModel++)
		{
			Triangle triangleModel = Triangle(molModel[indexModel],
											  molModel[indexModel + 1],
											  molModel[indexModel + 2]);

			//Calc the transformation to fit one triangle onto the other
			RigidTrans3 trans3 = triangleTarget->operator|(triangleModel);

			// match is a class that stores the correspondence list, eg.
			// pairs of atoms, one from eRigidTransach molecule, that are matching
			Match match;

			// apply rotation on each atom in the model molecule and
			// add the pairs of atoms (one from target and one from model)
			// that are close enough to the match list
			for (unsigned int i = 0; i < molModel.size(); i++)
			{
				Vector3 mol_atom = trans3 * molModel[i].position();

				// find close target molecule atoms using the hash
				HashResult<int> result;
				gHash.query(mol_atom, epsilon, result); // key is mol atom coordinate

				// check if the atoms in the result are inside the distance threshold
				// the hash is a cube shape, there can be atoms further that the threshold
				for (auto x = result.begin(); x != result.end(); x++)
				{
					float dist = mol_atom.dist(molTarget[*x].position());
					if (dist <= epsilon)
					{
						float score = (1 / (1 + dist));
						match.add(*x, i, score, score);
					}
				}
				result.clear();
			}

			//calculates transformation that is a little better than "rotation"
			match.calculateBestFit(molTarget, molModel);

			if (iMaxSize < match.size())
			{
				iMaxSize = match.size();
				rtransBest = match.rigidTrans();
				bestMatch = match;
			}
		}
	}
	//Export the alignment results
	std::cout << "Max Alignment Size: " << iMaxSize << std::endl;
	std::cout << "The achieved RMSD is: " << bestMatch.rmsd() << std::endl;
	std::cout << "Rigid Trans: " <<
			  RigidTrans3(Vector3(0, 0, 0), vectTargetMass) *
			  rtransBest *
			  RigidTrans3(Vector3(0, 0, 0), (-vectModelMass)) << std::endl;

	//Export the model molecule to PDB file
	std::ifstream file(argv[MODEL_INDEX]);
	Molecule<Atom> AtomsAfterTransformation = createUpdatedMolModel(&file, rtransBest);
	exportToPDB("transformed.pdb", AtomsAfterTransformation);

	//Calc the runtime of this program
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

}
