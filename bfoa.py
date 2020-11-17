# Bacterial Foraging Optimization Algorithm

import os, random, math, csv

class BFOA():

  def __init__(self, pop_size = 120, problem_size = 2, dimension = [-1, 1], elim_disp_steps = 1, repro_steps = 4, chem_steps = 40):
    self.step_index = 0
    self.run_id = os.urandom(6).encode('hex')

    # initializing variables received as arguments
    self.problem_size = problem_size
    self.dimension = [-1, 1]
    self.search_space = [self.dimension for x in range(self.problem_size)]
    
    self.pop_size = pop_size
    self.step_size = 0.1 # Ci

    self.elim_disp_steps = 1 # NED, configuring elimination-dispersal steps
    self.repro_steps = repro_steps # NRE, configuring reproduction steps

    self.chem_steps = chem_steps # NC, chemotaxis steps
    self.swim_length = 3 # NS, swim steps for cell
    self.p_eliminate = 0.25 # PED

    self.d_attr = 0.1 # attraction coefficients
    self.w_attr = 0.2
    self.h_rep = self.d_attr # coefficient of repulsion
    self.w_rep = 10

    # Setting up random generated population
    self.cells = [{'vector' : self.rnd_vector(self.search_space)} for x in range(self.pop_size)]
    print (self.search_space)

  # objective function
  def obj_function(self, vector):
    return sum(x**2.0 for x in vector)

  # function to create random number
  def rnd_vector(self, minmax):
    return [random.uniform(x[0], x[1]) for x in minmax]


  def generate_random_direction(self):
    return self.rnd_vector([self.dimension for x in range(self.problem_size)])


  def compute_cell_interaction(self, cell, d, w):
    '''
      Compare the present  cell to the other cells for attract or repel forces
    '''
    sum = 0.0

    for other_cell in self.cells:
      diff = 0.0
      for idx, i in enumerate( cell['vector'] ):
        diff += (cell['vector'][idx] - other_cell['vector'][idx])**2.0

      sum += d * math.exp(w * diff)
    return sum


  def attract_repel(self, cell):
    '''
      Compute the competing forces
    '''
    attract = self.compute_cell_interaction(cell, -self.d_attr, -self.w_attr)
    repel = self.compute_cell_interaction(cell, self.h_rep, -self.w_rep)
    return attract + repel


  def evaluate(self, cell):
    cell['cost'] = self.obj_function( cell['vector'] )
    cell['inter'] = self.attract_repel(cell)
    cell['fitness'] = cell['cost'] + cell['inter']
    return cell


  def tumble_cell(self, cell):
    step = self.generate_random_direction()

    vector = [None] * len(self.search_space)
    for idx, i in enumerate(vector):

      # Moving to next step in same direction as it is currently moving by step distance
      vector[idx] = cell['vector'][idx] + self.step_size * step[idx]

      # Checking for the condition if step takes out of environmnt and if so, then staying at the edge
      if vector[idx] < self.search_space[idx][0]: vector[idx] = self.search_space[idx][0]
      if vector[idx] > self.search_space[idx][1]: vector[idx] = self.search_space[idx][1]

    return {'vector' : vector}


  def save(self):
      #Adding columns for cell position in file
      with open('bfoa_'+ self.run_id +'.csv', 'a+eb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(["step_index", "x", "y", "cost", "inter", "fitness", "sum_nutrients"])

        for cell in self.cells:
          arr = [self.step_index, cell['vector'][0], cell['vector'][1], cell['cost'], cell['inter'], cell['fitness'], cell['sum_nutrients'] ]
          writer.writerow(arr)


  def chemotaxis(self):
    '''
      Best returns a cell instance
    '''
    best = None

    # Loop forchemotaxis steps
    for j in range(self.chem_steps):
      moved_cells = []

      # Looping through each cell present in the popolation 
      for cell_idx, cell in enumerate(self.cells):

        sum_nutrients = 0.0
        # Determine J of current cell position
        cell = self.evaluate(cell)

        # If this step give cell a lower energy result than the previous one
        if best is None or cell['cost'] < best['cost']: best = cell
        sum_nutrients += cell['fitness']

        # Loop for swim or tumble cell at some interval
        for m in range(self.swim_length):

          # Moving cell to new location
          new_cell = self.tumble_cell(cell)
          # Find J of the cell moved to new position 
          new_cell = self.evaluate(new_cell)

          # check if cell moved to new location has lower J value
          if cell['cost'] < best['cost']: best = cell
          # If not then try again
          if new_cell['fitness'] > cell['fitness']: break
          
          #  Save the total amount of food if new cell is better
          cell = new_cell
          sum_nutrients += cell['fitness']

        cell['sum_nutrients'] = sum_nutrients
        moved_cells.append( cell )

      print "  >> chemo=#{0}, f={1}, cost={2}".format(j, best['fitness'], best['cost'] )
      self.cells = moved_cells
      # capture these steps
      self.save()
      self.step_index += 1

    return best


  def search(self):
    '''
      Algorithm to iterate over a new random population
    '''
    best = None

    # Elimination dispersal step to delet cells and add new cell with low random probability
    for l in range(self.elim_disp_steps):

      # Reproduction Step: If a cell performed well over lifetime than it can contribute to next generation 
      for k in range(self.repro_steps):

        # Chemotaxis Step: cost of cells is estimated by the proximity to other cells and cells move along the manipulated cost surface one at a time
        # returns a single cell
        c_best = self.chemotaxis()

        # If the primary time, or if this reproduction step gave a lower energy cell
        if best is None or c_best['cost'] < best['cost']: best = c_best
        print " > best fitness={0}, cost={1}".format( best['fitness'], best['cost'] )

        # In reproduction, low health metric are eliminated and other healthy half population is doubled 
        self.cells = sorted(self.cells, key=lambda k: k['sum_nutrients'])
        lowest_cost_cells = self.cells[:self.pop_size/2]
        self.cells =  lowest_cost_cells + lowest_cost_cells

        # Save and increment index
        self.save()
        self.step_index += 1


      # Elimination-dispersal over each cell
      for cell in self.cells:
        if random.random() <= self.p_eliminate: cell['vector'] = self.rnd_vector(self.search_space)


      self.save()
      self.step_index += 1


    print "best :: ", best
    return best


if __name__ == "__main__":
  try:
    print "Starting BFOA ALGO.."
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    base_path = os.path.join(dname, 'bfoa_results')
    os.chdir(base_path)
  except RuntimeError:
    print "Create 'bfoa_results' folder in the directory"

  bfoa = BFOA(pop_size = 600, elim_disp_steps = 3, repro_steps = 4, chem_steps = 40, )
  best = bfoa.search()
  