'''
    Documentation of function parameters:
    @param <name> [<S.I. unit or type>] (<default value>): <Description>
'''

from math import pi


class Island:
    rho = 0.66
    L = 100

    def __init__(self, x, y, a, fixed=False, rho=None, L=None):
        '''
            @param x [nm]: X position of the island.
            @param y [nm]: Y position of the island.
            @param a [rad]: Angle by which the island is rotated.
            @param fixed [bool] (False): Whether the magnetization is fixed at <a>.
            @param rho [-] (Island.rho): The roundness of the island.
            @param L [nm] (Island.L): The major axis of the ellipses constituting the island.
        '''
        self.x = x
        self.y = y
        self.a = a
        self.fixed = fixed
        self.rho = rho
        self.L = L
    
    def init_geom(self):
        self.rho = Island.rho if self.rho is None else self.rho
        self.L = Island.L if self.L is None else self.L

    def set_geom(self, rho=None, L=None):
        self.rho = Island.rho if rho is None else rho
        self.L = Island.L if L is None else L
    
    def move(self, dx, dy):
        self.x += dx
        self.y += dy
    
    def __str__(self):
        return f"'x':{self.x}, 'y':{self.y}, 'a':{self.a}, 'fixed':{self.fixed}, 'rho':{self.rho}, 'L':{self.L}"


def generate_mumax3_inputfile(grid, islands, rho=None, L=None):
    '''
        @param grid [nm]: Size of a grid cell, in all dimensions.
        @param islands [list]: List of all Islands in the simulation.
        @param rho [-] (Island.rho): Roundness of all islands.
        @param L [nm] (Island.rho): Major axis of ellipses constituting the islands.

        If rho or L is specified in the Island __init__ call, then those are used. If not, then the values passed
            to generate_mumax3_inputfile are used. If those are not specified either, the default Island.rho and
            Island.L from the class Island are used.
    '''
    if rho is not None:
        for island in islands:
            if island.rho is None:
                island.rho = rho
    if L is not None:
        for island in islands:
            if island.L is None:
                island.L = L
    for island in islands:
        island.init_geom() # VERY IMPORTANT OTHERWISE rho AND L ARE None !!!
    
    #### Determine the simulation area
    # Simulation side lengths rounded to nearest power of two
    min_x = min_y = max_x = max_y = 0
    for island in islands:
        min_x = min(min_x, island.x-island.L/2)
        min_y = min(min_y, island.y-island.L/2)
        max_x = max(max_x, island.x+island.L/2)
        max_y = max(max_y, island.y+island.L/2)
    dx, dy = -(max_x + min_x)/2, -(max_y + min_y)/2 # How far to move the whole geometry to center it
    dx, dy = (dx//grid)*grid, (dy//grid)*grid # Make sure the position relative to the grid is kept the same
    for island in islands: island.move(dx, dy) # Move islands to be centered around (0,0)
    L_x = max(max_x+dx, -(min_x+dx))*2 # Since max_x+dx and min_x+dx should be centered, take the max of their abs to be certain
    L_y = max(max_y+dy, -(min_y+dy))*2
    Nx = 1
    while Nx < L_x/grid: Nx *= 2
    Ny = 1
    while Ny < L_y/grid: Ny *= 2

    #### Create the inputfile
    with open('many_islands_interaction.mx3template', 'r') as templateFile:
        text = ''.join([line for line in templateFile])
    # Grid
    text = text.replace(r'@{Nx}', str(int(Nx)))
    text = text.replace(r'@{Ny}', str(int(Ny)))
    text = text.replace(r'@{cellSize}', '%.2fe-9' % grid)
    # angle{n} (geometry angle) and a{n} (magnetization angle)
    angle_text = '\n'.join([('angle%d := %.10f' % (i+1, island.a)) for i, island in enumerate(islands)])
    text = text.replace(r'@{anglen}', angle_text)
    a_text = '\n'.join([('a%d := angle%d' % (i+1, i+1)) for i, island in enumerate(islands)])
    text = text.replace(r'@{an}', a_text)
    # islands and geometry
    islands_text = '\n'.join([('island%d := Ellipse(%.2fe-9, %.2fe-9)' % (i+1, island.L, island.L*island.rho)) for i, island in enumerate(islands)])
    islands_text += '\n' + '\n'.join([('island%d = island%d.Add(island%d.RotZ(Pi/2))' % (i+1,i+1,i+1)) for i, _ in enumerate(islands)])
    islands_text += '\n' + '\n'.join([('island%d = island%d.RotZ(angle%d).Transl(%.3fe-9, %.3fe-9, 0)' % (i+1,i+1,i+1, island.x, island.y)) for i, island in enumerate(islands)])
    text = text.replace(r'@{islandn}', islands_text)
    geometry_text = 'geometry := island1' + ''.join(['.Add(island%d)' % (i+2) for i, _ in enumerate(islands[1:])])
    text = text.replace(r'@{geometry}', geometry_text)
    # regions
    defregion_text = '\n'.join([('DefRegion(%d, island%d)' % (i+1, i+1)) for i, _ in enumerate(islands)])
    text = text.replace(r'@{DefRegion}', defregion_text)
    # external field for fixed islands
    extfield_text = '\n'.join([('B_ext.setRegion(%d, Vector(cos(angle%d), sin(angle%d), 0).Mul(fixationField))' % (i+1,i+1,i+1)) for i, island in enumerate(islands) if island.fixed])
    extfield_text += '\n' + '\n'.join([('m.setRegion(%d, Uniform(cos(angle%d), sin(angle%d), 0))' % (i+1,i+1,i+1)) for i, island in enumerate(islands) if island.fixed])
    text = text.replace(r'@{B_ext}', extfield_text)
    # table columns
    table_an_text = '\n'.join([('TableAddVar(a%d, "a%d", "rad")' % (i+1, i+1)) for i, island in enumerate(islands)])
    text = text.replace(r'@{TableAddan}', table_an_text)
    table_mregions_text = '\n'.join([('TableAdd(m.Region(%d))' % (i+1)) for i, _ in enumerate(islands)])
    text = text.replace(r'@{TableAddmRegions}', table_mregions_text)
    table_globals_text = ''
    if rho is not None:
        table_globals_text += 'TableAddVar(%.3f, "roundness", "")\n' % rho
    if L is not None:
        table_globals_text += 'TableAddVar(%.3f, "size", "m")\n' % L
    text = text.replace(r'@{TableAddGlobals}', table_globals_text)
    # loops
    loops_text = '\n'.join([('for a%d=angle%d; a%d < 2*Pi+angle%d; a%d += Pi/2 {' % (i+1,i+1,i+1,i+1,i+1)) for i, island in enumerate(islands) if not island.fixed])
    inside_text = '\n'.join([('    m.setRegion(%d, Uniform(1,0,0).rotZ(a%d))' % (i+1,i+1)) for i, island in enumerate(islands) if not island.fixed])
    text = text.replace(r'@{loops}', loops_text+'\n'+inside_text)
    text = text.replace(r'@{loops_closing_braces}', '}\n'*len([i for i in islands if not i.fixed]))
    
    with open('many_islands_interaction.mx3', 'w') as outFile:
        outFile.write(text)


if __name__ == "__main__":
    generate_mumax3_inputfile(2, [ # Basic half adder
        Island(-128, 0, 0),
        Island(0, 0, 0),
        Island(20, -180, pi/2, fixed=True)
    ], rho=0.66, L=100)
    # generate_mumax3_inputfile(2, [ # Basic half adder (with long wires that dont work)
    #     Island(-128, 0, 0),
    #     Island(0, 0, 0),
    #     Island(20, -156, pi/2, fixed=True),
    #     Island(-128, 128, 0),
    #     Island(90, 90, 0),
    #     Island(-128, 256, 0),
    #     Island(180, 180, 0)
    # ], rho=0.66, L=100)
    # generate_mumax3_inputfile(2, [ # Rhombic half adder shape (with smol wire)
    #     Island(0, 0, pi/4),
    #     Island(90, -90, 0),
    #     Island(-90, -90, 0),
    #     Island(0, -180, pi/4),
    #     Island(-180, -180, 0)
    # ], rho=0.66, L=100)
    # generate_mumax3_inputfile(2, [ # Small test for input on rhombic half adder
    #     Island(-90, -90, 0),
    #     Island(0, -180, pi/4),
    #     Island(-180, -180, 0)
    # ], rho=0.66, L=100)
