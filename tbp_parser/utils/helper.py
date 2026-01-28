import re


class Helper:
    """
    A collection of helper functions for parsing and analyzing mutations.
    These functions should not depend on any instance-specific data or require
    any import of other modules beyond standard libraries.
    """

    @staticmethod
    def get_position(mutation) -> list[int]:
        """This function recieves a mutation and returns the position as an integer
        Args:
            mutation (str): the mutation; can be either p.Met291Ile or p.Lys123_delArg125

        Returns:
            list[int]: the numerical position(s) of the mutation (in the above example, [291] or [123, 125])
        """
        pattern = r"-?\d+"
        match = re.findall(pattern, mutation)
        if len(match) > 0:
            return [int(x) for x in match]
        return [None]

    @staticmethod
    def get_mutation_genomic_positions(position, mutation) -> tuple[int, int]:
        """This function receives the genomic position and a mutation and returns the genomic position range as a list of integers

        Args:
            position(int): the genomic position of the mutation (e.g., 2000)
            mutation (str): the nucleotide mutation; can be either c.-33_327del or c.1693G>T

        Returns:
            tuple[int, int]: the start and stop genomic positions of the mutation (in the above example, [2000, 2360] or [2000, 2000])
        """
        pattern = r"-?\d+"
        match = re.findall(pattern, mutation)
        if len(match) == 1:
            return (position, position)
        elif len(match) == 2:
            return (position, position + (abs(int(match[0]) - int(match[1]))))
        return (None, None)

    @staticmethod
    def is_mutation_within_range(position, range_positions) -> bool:
        """Determines if a position is within a particular range

        Args:
            position (list[int]): either one or two positions
            range_positions (list[int] or list[list[int]]): the start and end regions of the range

        Returns:
            bool: true if the position is within the range_positions, false otherwise
        """
        try:
            if isinstance(range_positions[0], list):
                # check if the value is a list of lists; if so, check both lists
                return Helper.is_mutation_within_range(position, range_positions[0]) or Helper.is_mutation_within_range(position, range_positions[1])

            if len(position) > 1:
                # if the value is a list of two items, check if the position is within the range
                if any([x in range(range_positions[0], range_positions[1]) for x in position]):
                    return True
                if any([x in range(position[0], position[1]) for x in range_positions]):
                    return True

            # the position is a single item
            elif range_positions[0] <= position[0] <= range_positions[1]:
                return True

            return False
        except:
            return False