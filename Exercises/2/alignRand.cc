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


#define NUMBER_OF_ARGS 4
#define EPSILON_INDEX 1
#define TARGET_INDEX 2
#define MODEL_INDEX 3


void moveToCenter(Molecule<Atom> *molModel)
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
}

Molecule<Atom> createMolModel(std::ifstream *sourceFile, RigidTrans3 rtrans)
{
	Molecule<Atom> molModelForExport;
	molModelForExport.readPDBfile(*sourceFile, PDB::AllSelector());

	// calculate center of mass
	moveToCenter(&molModelForExport);

	//Preform rotations
	Match match;
	for (unsigned int i = 0; i < molModelForExport.size(); i++)
	{
		Vector3 mol_atom = rtrans * molModelForExport[i].position(); // rotate
		molModelForExport[i].update(mol_atom);
	}

	return molModelForExport;
}

void exportToPDB(std::string fileName, Molecule<Atom> molModelToExport)
{
	std::ofstream of(fileName);
	if (of.is_open())
	{
		of << molModelToExport << std::endl;
		of.flush();
		of.close();
		std::cout << "wrote the file successfully!" << std::endl;
	}
	else
	{
		std::cerr << "Failed to open file transformed.pdb" << std::endl;
		exit(-1);
	}
}

bool checkIfRNA(std::string fileName)
{
	std::string RNA = " RNA ";
	std::ifstream file(fileName);
	std::string line;
	getline(file, line);
	return line.find(RNA, 0) != std::string::npos;
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

	std::cout << "Given epsilon: " << epsilon << std::endl;


	// read the two files into Molecule
	Molecule<Atom> molModel, molTarget;

	std::ifstream fileTarget(argv[TARGET_INDEX]);
	if (!fileTarget)
	{
		std::cout << "File " << argv[TARGET_INDEX] << " does not exist." << std::endl;
		return 0;
	}
	else
	{
		std::cout << "File " << argv[TARGET_INDEX] << " was successfully found." << std::endl;
	}

	std::ifstream fileModel(argv[MODEL_INDEX]);
	if (!fileModel)
	{
		std::cout << "File " << argv[MODEL_INDEX] << " does not exist." << std::endl;
		return 0;
	}
	else
	{
		std::cout << "File " << argv[MODEL_INDEX] << " was successfully found." << std::endl;
	}

	if (checkIfRNA(argv[MODEL_INDEX]))
	{
		molTarget.readPDBfile(fileTarget, PDB::PSelector());
		molModel.readPDBfile(fileModel, PDB::PSelector());
	}
	else
	{
		molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
		molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
	}

	// calculate center of mass
	Vector3 vectModelMass(0, 0, 0);
	for (unsigned int i = 0; i < molModel.size(); i++)
	{
		vectModelMass += molModel[i].position();
	}
	vectModelMass /= molModel.size();

	Vector3 vectTargetMass(0, 0, 0);
	for (unsigned int i = 0; i < molTarget.size(); i++)
	{
		vectTargetMass += molTarget[i].position();
	}
	vectTargetMass /= molTarget.size();

	// transform the molecules to the center of the coordinate system
	molModel += (-vectModelMass);
	molTarget += (-vectTargetMass);

	if (molModel.empty() || molTarget.empty())
	{
		std::cout << "No molecules were found in PDB file" << std::endl;
		return -1;
	}

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

	// now we try random rotations and choose the best alignment from random rotations
	// פה זה החלק שבו אנחנו צריכים להתערב בקוד
	unsigned int iMaxSize = 0;
	Match bestMatch;
	RigidTrans3 rtransBest;

	unsigned int transCounter = 0;
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
//			Matrix3 rotation = trans3.rotation(); //get the rotation matrix
			transCounter++;


			// match is a class that stores the correspondence list, eg.
			// pairs of atoms, one from eRigidTransach molecule, that are matching
			Match match;

			// apply rotation on each atom in the model molecule and
			// add the pairs of atoms (one from target and one from model)
			// that are close enough to the match list
			for (unsigned int i = 0; i < molModel.size(); i++)
			{
				Vector3 mol_atom = trans3 * molModel[i].position();
//				Vector3 mol_atom = rotation * molModel[i].position(); // rotate

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

	std::cout << "Max Alignment Size: " << iMaxSize << std::endl;
	std::cout << "The achieved RMSD is: " << bestMatch.rmsd() << std::endl;
	std::cout << "Rigid Trans: " <<
			  RigidTrans3(Vector3(0, 0, 0), vectTargetMass) *
			  rtransBest *
			  RigidTrans3(Vector3(0, 0, 0), (-vectModelMass)) << std::endl;

	std::ifstream file(argv[MODEL_INDEX]);

	Molecule<Atom> AtomsAfterTransformation = createMolModel(&file, rtransBest);
	exportToPDB("transformed.pdb", AtomsAfterTransformation);

	std::ifstream file2(argv[TARGET_INDEX]);
	Molecule<Atom> molModelForExport;
	molModelForExport.readPDBfile(file2, PDB::AllSelector());
	moveToCenter(&molModelForExport);
	exportToPDB("transformed2.pdb", molModelForExport);

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
	std::cout << "# of transformations: " << transCounter << std::endl;

}
