#################################################################
# FILE : nonogram.py
# WRITER : Gold Nimrod , nimrod
#          Eldar Michal , michalel
# EXERCISE : intro2cs2 ex8 2021
# DESCRIPTION: The code for ex8 to solve a nonogram puzzle
# NOTES: We chose to return an unknown cell in the intersection rows, if there
#        on a certain column we have unknown cells and a specific color.
#        We chose this option because we felt that this implementation is more
#        efficient than the other way around (if we found at least one unknown
#        cell, we exit immediately).
#################################################################

# Color consts
BLACK = 1
WHITE = 0
UNKNOWN = -1

# Constraints indexes
ROW = 0
COLUMN = 1


def get_first_block_sequence(index, blocks):
    """
    Creates a sequence representing the first block in the blocks constraints
    :param index: The end to end the block from
    :param blocks: The list of block constraints
    :return: A list representing the sequence of the first block
    """
    block = index * [WHITE] + blocks[0] * [BLACK]
    # If it isn't the last block - add a white separator
    if len(blocks) > 1:
        block += [WHITE]
    return block


def constraint_satisfactions(n, blocks):
    """
    Returns all the valid rows of a specific length which satisfy the block
    constrains
    :param n: The length of the row
    :param blocks: The list of block constraints
    :return: A list with all the valid row options
    """
    # If there are no block constrains - return a row with only white cells
    if len(blocks) == 0:
        return [n * [WHITE]]

    options = []
    # Run on all the optional starting cells of the first block
    for i in range(n):
        if i + blocks[0] > n:
            return options

        start_block = get_first_block_sequence(i, blocks)
        # Get all the options of the rest of the row
        for option in constraint_satisfactions(n - len(start_block),
                                               blocks[1:]):
            options.append(start_block + option)
    return options


def get_block_start_option(row, start):
    """
    Returns the first last_guess_column that can be the end of a block -
    Either a black cell or an unknown one
    :param row: The row to search in
    :param start: The last_guess_column to end the search from
    :return: The first last_guess_column that is black or unknown.
             None if there are only white cells left
    """
    for i in range(start, len(row)):
        if row[i] in [BLACK, UNKNOWN]:
            return i
    return None


def empty_unknowns_in_range(row, start, end):
    """
    Paint in white all the unknowns in a certain range - Changes the given row
    :param row: The row to change
    :param start: The starting last_guess_column (inclusive)
    :param end: The ending last_guess_column (exclusive)
    :return: None
    """
    for i in range(start, end):
        if row[i] == UNKNOWN:
            row[i] = WHITE


def fill_end_with_white(row, start):
    """
    Fills squares with white from start index to the end
    (means no more blocks left to fill)
    :param row: row to fill (list)
    :param start: index in the row (int)
    :return: row, completed with white.
    """
    for i in range(start, len(row)):
        if row[i] == BLACK:
            return []
        if row[i] == UNKNOWN:
            row[i] = WHITE
    return row


def add_block_separator(row, separator_index):
    """
    Adds a separator (if needed) after the block
    :param row: The row to change
    :param separator_index: The last_guess_column of the separator to paint
    :return: None
    """
    if separator_index < len(row):
        if row[separator_index] == UNKNOWN:
            row[separator_index] = WHITE


def is_valid_block_end(row, end_block_index):
    """
    Checks if the end of the block is valid - meaning it isn't colored
    :param row: The row to check
    :param end_block_index: The end of the block to check
    :return: True if the end is valid. False otherwise
    """
    # If the ending is the last square in the row, we are okay
    if end_block_index < len(row) - 1:
        if row[end_block_index] == BLACK:
            return False
    return True


def fill_block(row, block, start):
    """
    Fills a block in a row, at the given index, if possible
    :param row: The row to change
    :param block: The block we want to add
    :param start: The wanted start position of the block
    :return: The row if the block was added.
             An empty list if the board can't be added in the given position
    """
    separator_index = start + block
    # If the last cell in the block is white, or if after the block we have
    # a black cell, we can't put the block here
    if row[separator_index - 1] == WHITE or \
            not is_valid_block_end(row, separator_index):
        return []

    for i in range(start, separator_index):
        # If we reached a white cell, this block isn't valid
        if row[i] == WHITE:
            return []
        if row[i] == UNKNOWN:
            row[i] = BLACK
    # Adds a white cell separator after the block (if needed)
    add_block_separator(row, separator_index)
    return row


def row_variations_helper(row, blocks, current_index):
    """
    Gets row, and blocks to fill. Returns all the valid variations (options) of
    filling the raw with given blocks.
    :param row: row to fill (list)
    :param blocks: blocks to fill in row (list of lists)
    :param current_index: index to start fill from.
    :return: list or row variations.
    """
    if not row:
        return [[]]
    if not blocks:
        filled_row = row[:]
        filled_row = fill_end_with_white(filled_row, current_index)
        return [filled_row] if filled_row else []

    variations = []
    filled_row = row[:]
    continue_searching = True
    start_of_block = get_block_start_option(filled_row, current_index)

    while continue_searching and start_of_block is not None:
        # If we don't have enough space to fill the block
        if blocks[0] > len(filled_row) - start_of_block:
            return variations
        empty_unknowns_in_range(filled_row, 0, start_of_block)

        # If The first valid end is a black cell, our block can't be after it
        if filled_row[start_of_block] == BLACK:
            continue_searching = False

        filled_row = fill_block(filled_row, blocks[0], start_of_block)
        # Only if we managed to fill the block
        if filled_row:
            block_end = start_of_block + blocks[0]
            variations += row_variations_helper(filled_row, blocks[1:],
                                                block_end)
        filled_row = row[:]
        start_of_block = get_block_start_option(filled_row, start_of_block + 1)
    return variations


def row_variations(row, blocks):
    """
    Check all row variations for filling given blocks in row.
    using row_variations_helper function
    :param row: row to fill
    :param blocks: blocks to fill in row
    :return: list of variations.
    """
    return row_variations_helper(row, blocks, 0)


def check_col(rows, col_index):
    """
    checks if specific column has the same value in all squares.
    If it has - returns the value. else, returns -1
    :param rows: rows to check
    :param col_index: index of column
    :return: common value
    """
    for i in range(1, len(rows)):
        if rows[i][col_index] != rows[i - 1][col_index]:
            return UNKNOWN
    return rows[0][col_index]


def intersection_row(rows):
    """
    returns an intersection of rows. If all squares value in the same index
    (column index) are equals - fills the intersection row with the value.
    else - fill with -1.
    :param rows: rows for intersection
    :return: intersection row
    """
    if len(rows) <= 1:
        return rows[0]

    rows_intersection = []
    for i in range(len(rows[0])):
        rows_intersection.append(check_col(rows, i))
    return rows_intersection


def rotate_board(board):
    """
    Rotates a board
    :param board: The board to rotate
    :return: the rotated board
    """
    rotated_board = []
    row_length = len(board[0])
    col_length = len(board)

    for i in range(row_length):
        row = []
        for j in range(col_length):
            row.append(board[j][i])
        rotated_board.append(row)
    return rotated_board


def solve_easy_nonogram(constraints):
    """
    solves nonogram only from given constraints, without making guesses.
    :param constraints: constraints for board
    :return: board's solution
    """
    board = create_board(constraints)
    board = solve_board(constraints, board)

    return board


def solve_board(constraints, board):
    """
    Solves board according to given constraints. Runs over each row and solve
    it using solve_by_rows functions. Ends when no more conclusions from
    constraint left to conclude.
    :param constraints: board constraints
    :param board: list of board's rows
    :return: solved board if board can be solved. else None.
    """
    if not board:
        return None

    row_constraints = constraints[ROW]
    col_constraints = constraints[COLUMN]

    changed = True
    while changed:
        row_solved = solve_by_rows(board, row_constraints)
        if not row_solved:
            return None
        board, row_changed = row_solved
        board = rotate_board(board)
        col_solved = solve_by_rows(board, col_constraints)
        if not col_solved:
            return None
        board, col_changed = col_solved
        changed = row_changed or col_changed
        board = rotate_board(board)
    return board


def solve_by_rows(board, row_constraints):
    """
    Tries solving a board according the the row constraints, without guessing
    :param board: The board to solve
    :param row_constraints: The constraints of the rows
    :return: The most solved state of the board.
             None if it found out that the board couldn't be solved
    """
    changed = False
    # Get the options of each row in the board
    for i in range(len(row_constraints)):
        options = row_variations(board[i], row_constraints[i])

        # If there is no solution for the board, the board has no solution too
        if not options:
            return None

        new_row = intersection_row(options)
        # If we received new information, continue solving
        if new_row != board[i]:
            changed = True
        board[i] = new_row
    return board, changed


def create_board(constraints):
    """
    Creates a starting board according to the constraints, while filling the
    intersections of all the rows.
    :param constraints: The constraints of the board to create
    :return: The initialized board with the intersections of the rows
    """
    board = []

    for constraint in constraints[ROW]:
        row_options = constraint_satisfactions(len(constraints[COLUMN]),
                                               constraint)
        board.append(intersection_row(row_options))
    return board


def find_first_unknown(board, row_index, col_index):
    """
    Finds the first occurrence of an unknown cell
    :param board: The board to search in
    :param row_index: The index of the row to start searching from
    :param col_index: The index of the column to start searching from
    :return: A tuple consisting the position of the first unknown cell found
             None if there are no unknowns in the board
    """
    # In the first row - search just from the row index
    for j in range(col_index, len(board[0])):
        if board[row_index][j] == UNKNOWN:
            return row_index, j
    # In the next rows - search from the start of each row
    for i in range(row_index + 1, len(board)):
        for j in range(len(board[0])):
            if board[i][j] == UNKNOWN:
                return i, j
    # If we didn't find any unknown cells
    return None


def solve_nonogram_helper(constraints, board, last_guess_row,
                          last_guess_column, solutions):
    """
    Solves a specific nonogram board, while guessing if needed
    :param constraints: The list of the board's constraints
    :param board: The board to solve
    :param last_guess_row: The row of the last cell we guessed
    :param last_guess_column: The column of the last cell we guessed
    :param solutions: The list of the board's solutions - changes it inside
    :return: The list of board solutions
    """
    solved_board = solve_board(constraints, board)
    if not solved_board:
        return solutions

    first_unknown_position = find_first_unknown(solved_board, last_guess_row,
                                                last_guess_column)
    # If we couldn't solve the board entirely, we need to guess
    if first_unknown_position:
        row = first_unknown_position[ROW]
        col = first_unknown_position[COLUMN]
        # Guess both white and black for the cell,
        # Then try solving the new boards
        for color in [WHITE, BLACK]:
            board_guess = solved_board[:]
            board_guess[row][col] = color
            solve_nonogram_helper(constraints, board_guess, row, col,
                                  solutions)
    # If we reached a solution - add it
    else:
        solutions.append(solved_board)
    return solutions


def solve_nonogram(constraints):
    """
    Solves a nonogram while guessing if needed
    :param constraints: The list of the board's constraints
    :return: A list containing all the board's solutions
    """
    board = create_board(constraints)
    solutions = solve_nonogram_helper(constraints, board, 0, 0, [])

    return solutions


if __name__ == "__main__":
    n_10_15 = [[[4], [3, 1], [2, 1, 1], [1, 2, 1], [2, 2, 1, 1], [3, 1, 1, 1], [1, 1, 1, 2], [3, 5],
                [1, 1, 2], [2, 2, 1], [1, 3, 1], [3, 1], [1, 1], [4], [4]],
               [[5, 2, 2], [3, 2, 1, 1, 3], [2, 1, 5, 2], [1, 1, 3, 1, 2], [1, 2, 8], [2, 1, 1, 1, 1], [1, 3],
                [1, 1, 4], [1, 2], [2]]]

    n_20_20 = [
        [[7, 7], [8, 7], [11, 3], [4, 1, 3], [5, 1, 2, 6], [3, 4, 3], [5, 3, 3, 2], [4, 3, 6], [9, 2, 6], [11, 1, 2, 1],
         [1, 5, 4, 2], [1, 4, 2], [11], [4, 7], [4, 4], [3, 1], [3], [3], [3], [3]],
        [[1, 3], [3, 4], [5, 5], [5, 5, 6], [11, 7], [7, 3, 8], [3, 2, 7], [3, 1, 2, 2], [1, 3, 3, 1], [3, 3, 5],
         [3, 2, 1, 5], [3, 5, 5], [2, 1, 2, 4], [5, 1, 4], [5, 8], [5, 7], [1, 2, 2], [2, 2], [5], [4]]]

    n_25_30 = [
        [[3], [4], [4], [5], [10], [13], [1, 15], [17, 2], [16, 4], [15, 6], [3, 10, 2, 3], [2, 10, 2, 4], [2, 6, 5, 5],
         [2, 5, 2, 7], [2, 4, 10, 1], [2, 3, 1, 2, 6, 1], [2, 3, 2, 2, 4, 1], [3, 2, 4, 2, 4, 1], [6, 1, 2, 1, 3, 1],
         [5, 3, 1, 3, 1], [4, 1, 2, 3, 1], [1, 1, 2, 2, 1], [1, 1, 3, 1, 1], [3, 3, 6], [2, 2, 8], [2, 2, 10],
         [2, 1, 2, 5], [3, 1, 4], [3], [2]],
        [[7], [15], [5, 5], [4, 3, 2], [3, 10, 5], [14, 5], [12, 1, 1], [10, 2, 4], [10, 2, 1, 5], [9, 6, 2], [8, 3, 2],
         [9, 2], [12, 5], [12, 4], [11, 4], [10, 3, 3], [1, 5, 2, 1, 3], [4, 3, 2, 3], [2, 3, 5, 4], [4, 10, 3],
         [15, 4], [13, 6], [7, 1, 6], [3, 2, 6], [6, 5]]]
    import time

    start_time = time.time()
    solve_nonogram(n_20_20)
    print("--- %s seconds ---" % (time.time() - start_time))

