from game import Agent
import api
import itertools
from pacman import Directions

''' 
The following Grid class has been used as-is courtesy of Simon Parsons with nothing added but some
unused functionality removed. This class and all other functions borrowed from Simon Parsons are only used
for printing out utilities to the console, my MDP solver uses dictionaries rather than this grid object.
'''

# A class that creates a grid that can be used as a map
#
# The map itself is implemented as a nested list, and the interface
# allows it to be accessed by specifying x, y locations.
class Grid:

    # Constructor
    #
    # Note that it creates variables:
    #
    # grid:   an array that has one position for each element in the grid.
    # width:  the width of the grid
    # height: the height of the grid
    #
    # Grid elements are not restricted, so you can place whatever you
    # like at each location. You just have to be careful how you
    # handle the elements when you use them.
    def __init__(self, width, height):
        self.width = width
        self.height = height
        subgrid = []
        for i in range(self.height):
            row = []
            for j in range(self.width):
                row.append(0)
            subgrid.append(row)

        self.grid = subgrid

    # The display function prints the grid out upside down. This
    # prints the grid out so that it matches the view we see when we
    # look at Pacman.
    def prettyDisplay(self):
        for i in range(self.height):
            for j in range(self.width):
                # print grid elements with no newline
                print self.grid[self.height - (i + 1)][j],
            # A new line after each line of the grid
            print
            # A line after the grid
        print

    # Set and get the values of specific elements in the grid.
    # Here x and y are indices.
    def setValue(self, x, y, value):
        self.grid[y][x] = value

    def getValue(self, x, y):
        return self.grid[y][x]

    # Return width and height to support functions that manipulate the
    # values stored in the grid.
    def getHeight(self):
        return self.height

    def getWidth(self):
        return self.width


"""
An agent class that plays Pacman in a non-deterministic environment using an MDP solver based on value iteration.

The agent should beat the benchmarks with its current configuration, however if it fails please vary the 'r' argument 
of ghostRadius() on line 186 to '3' for mediumClassic or '2' for smallGrid. 

"""

class MDPAgent(Agent):

    def __init__(self):
        print "Here we go again..."
        # initialise dictionary variables for use in value iteration
        self.rewardDict = {}
        self.stateDict = {}
        self.utils = {}

        # variables for use in assigning rewards to map locations. Here and elsewhere long iterables such as walls and
        # other map features are stored as sets for speedy iteration
        self.walls = set()
        self.grid = set()

    def registerInitialState(self, state):
        """ allows initial stock-taking of environment at start-up via 'state' argument """

        # determine a dictionary of possible states for every walkable square in the map
        self.stateMapping(state)

        # create a Grid instance for printing out utilities
        self.makeMap(state)  # courtesy of mapAgents.py via Simon Parsons

        self.addWallsToMap(state)  # courtesy of mapAgents.py via Simon Parsons

        self.walls = set(api.walls(state))

        # determine (x, y) co-ordinate pairs for the map grid and store in a class attribute
        # for later use in mapping out actions, states and rewards within the game
        corners = api.corners(state)
        BL = (0, 0)
        BR = corners[1]
        TL = corners[2]
        map_width = range(BL[0], BR[0])
        map_height = range(BL[1], TL[1])
        self.grid = set((x, y) for x in map_width for y in map_height)

    ''' Courtesy of mapAgents.py via Simon Parsons '''
    # Make a map by creating a grid of the right size
    def makeMap(self, state):
        corners = api.corners(state)
        height = self.getLayoutHeight(corners)
        width = self.getLayoutWidth(corners)
        self.map = Grid(width, height)

    ''' Courtesy of mapAgents.py via Simon Parsons '''
    # Functions to get the height and the width of the grid.
    #
    # We add one to the value returned by corners to switch from the
    # index (returned by corners) to the size of the grid (that damn
    # "start counting at zero" thing again).
    def getLayoutHeight(self, corners):
        height = -1
        for i in range(len(corners)):
            if corners[i][1] > height:
                height = corners[i][1]
        return height + 1

    def getLayoutWidth(self, corners):
        width = -1
        for i in range(len(corners)):
            if corners[i][0] > width:
                width = corners[i][0]
        return width + 1

    ''' Courtesy of mapAgents.py via Simon Parsons '''
    # Functions to manipulate the map.
    #
    # Put every element in the list of wall elements into the map
    def addWallsToMap(self, state):
        walls = api.walls(state)
        for i in range(len(walls)):
            self.map.setValue(walls[i][0], walls[i][1], '|    |')


    def rewardMapping(self, state):
        """ Populate a dictionary of reward values for each traversable square in the map """

        walls = self.walls
        food = set(api.food(state))
        ghostStates = api.ghostStates(state)
        capsules = set(api.capsules(state))

        # initialise all states to reward of -1, the default value for a blank square
        self.rewardDict = {key: -1 for key in self.grid if key not in walls}

        # initialise all states' utilities to 0 whilst we're at it
        self.utils = {key: 0 for key in self.grid if key not in walls}

        # update the reward dictionary with reward values at locations of food and capsules
        foodDict = {k: 10 for k, v in self.rewardDict.items() if k in food}
        self.rewardDict.update(foodDict)
        capsuleDict = {k: 20 for k, v in self.rewardDict.items() if k in capsules}
        self.rewardDict.update(capsuleDict)

        # loop through the ghost locations and their statuses indicating whether they're scared; if they
        # are scared, ignore them. Otherwise, assign negative rewards to both: a) their locations and b)
        # a radius of variable width around them.
        for j in ghostStates:
            if j[0] in self.rewardDict.keys():
                if j[1] == 0:
                    self.rewardDict[j[0]] = -50

                    ''' adjust size of exclusion zone around ghosts with 3rd argument below '''
                    ghostNeighbours = self.ghostRadius(state, j[0], 2)

                    ghostRadius = {k: -25 for k, v in self.rewardDict.items() if k in ghostNeighbours}
                    self.rewardDict.update(ghostRadius)

    def stateMapping(self, state):
        """
        create a dictionary that maps all walkable grid squares to all of the states it's possible to reach from
        them based on the set of available actions at that square
        """
        walls = set(api.walls(state))

        stateDict = dict.fromkeys(self.rewardDict.keys())

        # loop through all squares in the map, for each one assign the possible states that could result from
        # carrying out the actions 'North', 'South', 'East' or 'West' at that point under the stochastic transition
        # model. Assign possible state values so that the intended result of each action is first.
        for i in stateDict.keys():
            tmp = self.neighbours(i)
            stateDict[i] = {'North': [tmp[3], tmp[0], tmp[2]],
                            'South': [tmp[1], tmp[0], tmp[2]],
                            'East': [tmp[0], tmp[3], tmp[1]],
                            'West': [tmp[2], tmp[3], tmp[1]],
                            }
            # if any of the possible states for a grid square represent walls, overwrite them with the initial grid
            # square value instead, i.e. if trying to move into a wall remain in position.
            for a, b in stateDict[i].items():
                for s in b:
                    if s in walls:
                        b[b.index(s)] = i

        self.stateDict = stateDict

    def ghostRadius(self, state, ghosts, r, next=None):
        """ Recursively creates a radius equal to 'r' moves around ghosts. Returns a set of (x, y) co-ords. """
        ghostLocs = api.ghosts(state)
        walls = set(api.walls(state))

        if r < 1:
            raise ValueError("'r' value for ghostRadius must be greater than 1")

        # generate neighbours of ghost locations on first call of function
        if next is None:
            next = []
            ghostNeighbours = self.neighbours(ghosts)
            ghostNeighbours = [i for i in ghostNeighbours if i not in walls]
            ghostNeighbours = [i for i in ghostNeighbours if i not in ghostLocs]

        # generate neighbours from the results of the function's previous pass
        if next:
            ghostNeighbours = [self.neighbours(i) for i in ghosts]
            ghostNeighbours = itertools.chain.from_iterable(ghostNeighbours)
            ghostNeighbours = [i for i in ghostNeighbours if i not in walls]
            ghostNeighbours = [i for i in ghostNeighbours if i not in ghostLocs]

        # return a final set of ghost neighbour locations if on the last pass
        if r == 1:
            next.append(set(ghostNeighbours))
            final = [list(i) for i in next]
            final = set(itertools.chain.from_iterable(final))
            return final

        # else decrement 'r' by one and call the function on itself
        else:
            r = r-1
            next.append(set(ghostNeighbours))
            return self.ghostRadius(state, set(ghostNeighbours), r, next)

    def valueIteration(self, state):
        """ implementation of the Value Iteration algorithm for solving an MDP """

        # empirically determined optimal gamma value
        gamma = 0.6

        # error bound for termination condition
        epsilon = 0.001

        # inputs representing various components of the Bellman equation
        states = self.stateDict  # set of states, S
        gridVals = self.rewardDict  # set of rewards, R(s)
        utilDict = self.utils  # set of utilities initialised at 0, U

        # loop indefinitely until convergence as determined by break condition
        while True:
            # delta is the change in utility for any one state from one iteration to the next; reset to 0 at start
            # of each iteration.
            delta = 0

            # for each square in the grid, register the current utility value for later assessment and generate a
            # temporary dictionary structure to store new utilities.
            for square, utility in utilDict.items():
                U = utility
                tmp_utils = {}

                # assign a new utility value to each grid square based on the reward at that square and the utilities of
                # the set of possible states from that square
                for direction, state in states[square].items():

                    U_s = gridVals[square] + gamma * (
                                0.8 * utilDict[state[0]] + 0.1 * utilDict[state[1]] + 0.1 * utilDict[state[2]])
                    tmp_utils[direction] = U_s

                # take the maximum expected utility value for the set of actions A(s) at each square, assign it to the
                # utility dictionary, U.
                utilDict[square] = max(tmp_utils.values())

                # find the difference in utility values between any state between iterations, check whether it satisfies
                # the termination condition and return the dictionary of final utilities if so
                delta = max(delta, abs(utilDict[square] - U))
            if delta < epsilon * (1 - gamma) / gamma:
                return utilDict

    def neighbours(self, id):
        """ returns a list of one-move neighbour locations of the 'id' argument """
        (x, y) = id
        results = [(x + 1, y), (x, y - 1), (x - 1, y), (x, y + 1)]  # E, S, W, N
        return results

    def bestMove(self, state, gridVals):
        """
        Based on Pacman's current position, return the neighbouring co-ordinate location associated with the
        maximum expected utility
        """
        walls = set(api.walls(state))
        loc = api.whereAmI(state)
        # generate neighbours of Pacman's location
        possibleStates = [i for i in self.neighbours(loc) if i not in walls]
        # the utilities of those neighbours
        U_states = [gridVals[i] for i in possibleStates]
        # the location of the maximum utility
        bestMove = possibleStates[U_states.index(max(U_states))]
        return bestMove

    def singleMove(self, location, neighbour):
        """ Convert (x, y) location of maximum expected utility neighbour from bestMove into a direction for
        getAction()
        """
        bearing = tuple(x - y for x, y in zip(location, neighbour))
        if bearing == (-1, 0):
            return Directions.EAST
        if bearing == (1, 0):
            return Directions.WEST
        if bearing == (0, -1):
            return Directions.NORTH
        if bearing == (0, 1):
            return Directions.SOUTH

    def getAction(self, state):

        start = api.whereAmI(state)

        # generate rewards for map locations in each game cycle to reflect changing conditions
        self.rewardMapping(state)

        # create the dictionary of states at the start of the game
        if not self.stateDict:
            self.stateMapping(state)

        # get converged utilities from value iteration
        gridVals = self.valueIteration(state)

        # get and assign the utility values and their associated locations to the Grid object, then print this map
        # representation to the console
        valKeys = gridVals.keys()
        valVals = gridVals.values()

        for i in range(len(gridVals)):
            y = valKeys[i][0]
            x = valKeys[i][1]
            val = '{:3.2f}'.format(valVals[i])
            self.map.setValue(y, x, val)

        self.map.prettyDisplay()

        # optimal policy for Pacman: get the location of the optimal utility from his current position
        sPrime = self.bestMove(state, gridVals)

        # required argument for makeMove()
        legal = api.legalActions(state)

        # return a move based on the direction returned by singleMove() and sPrime
        return api.makeMove(self.singleMove(start, sPrime), legal)

