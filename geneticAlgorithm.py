#!/usr/bin/env python3

# Please feel free to modify any part of the code or add any functions, or modify the argument
# of the given functions. But please keep the name of the given functions

# Please feel free to import any mathematical libraries you need.

# You are required to finish the genetic_algorithm function, and you may need to complete crossover, mutate and select.

import matplotlib.pyplot as plt
import random
from PIL import Image

def get_children(parents, probability_crossover, n_children=1):
    '''
    get a list of children based on crossover probability and number of children per parent pair
    '''
    shuffled = random.sample(parents, len(parents))
    children = []
    for i in range(0, len(shuffled), 2):
        for _ in range(n_children):
            children.append(copulate(probability_crossover, shuffled[i], shuffled[i+1]))
    return children

def copulate(probability_crossover, *parents):
    '''
    generate a child from multiple parents
    '''
    parent_index = random.randrange(len(parents))
    child = ''
    for i in range(30):
        child += parents[parent_index][i]
        if random.random() < probability_crossover:
            parent_index = (parent_index + 1) % len(parents)
    return child

def mutate(ants, probability_mutation):
    '''
    ants : list[str] -- ants to mutate
    probability_mutation : number -- how likely each character is to randomly mutate

    '''
    new_ants = []
    mutations = 0
    for ant in ants:
        new_ant = ''
        for i, c in enumerate(ant):
            if random.random() < probability_mutation:
                if i % 3 == 0:
                    new_ant += str(random.randint(1, 4))
                else:
                    new_ant += str(random.randint(0, 9))
            else:
                new_ant += c
        if ant != new_ant:
            mutations += 1
        new_ants.append(new_ant)
    return new_ants, mutations

def select_parents(old_gen, n_parents, auto_top_n=-1):
    '''
    old_gen : list[str] -- population to select from
    n_parents : int -- how many parents to select

    return : list[str] -- population of ants selected
    '''
    if auto_top_n == -1:
        auto_top_n = n_parents//4

    #preselect the top n performers (default is a quarter of the parents to be selected)
    chosen_indices = set(range(auto_top_n))
    indices = list(range(len(old_gen)))

    #lottery for the rest
    #higher weights for higher fitness
    weights = [(population_size - i)**2  for i in indices]
    for index in range(auto_top_n):
        weights[index] = 0

    remaining = n_parents
    while len(chosen_indices) != n_parents:
        remaining = n_parents - len(chosen_indices)
        for selected_index in random.choices(indices, weights, k=remaining):
            if weights[selected_index]:
                weights[selected_index] = 0
                chosen_indices.add(selected_index)
    # print(chosen_indices

    return [old_gen[i] for i in chosen_indices]

def select_elite(ants, fitness, n_elites):
    '''
    return n elite ants from a population. 
    fitness is included as a parameter for future extensibility.
    '''
    return ants[0:n_elites]

def genetic_algorithm(population_size, food_map, map_size, max_generation, probability_crossover, probability_mutation, 
        parent_proportion=0.6, n_children=3, greedy_proportion=0.1):
    analytics = {}
    analytics['best_fitness_per_gen'] = []


    # make sure parameters are valid
    n_parents = round(parent_proportion*population_size)
    n_elites = population_size - n_children * n_parents // 2
    try:
        assert (n_children*n_parents//2 + n_elites) == population_size
    except AssertionError:
        print('Error: n_children*n_parents//2 must be less than the population size')
        exit(1)


    # initialize and generate heatmaps for population
    print('==INITIAL POPULATION (generation 0)==')
    print('generating population...')
    ants = initialize_population(population_size, greedy_proportion)
    # print('index -- genes -- fitness')
    fitness = {}
    for i, ant in enumerate(ants):
        trial, fitness[ant], _ = ant_simulator(food_map, map_size, ant, heatmap=False)
        # display_trials(trial, f"initial/initial-{i}.txt")
        # display_trial_img(heatmap, map_size, f'initial/initial-{i}.png')
        # print(i, ant, fitness[ant])
    # print('heatmaps of initial population saved to /initial')
    print()


    ants = sorted(ants, key=lambda x: fitness[x], reverse=True)

    # run generations
    for generation_n in range(1, max_generation + 1):
        print(f'==GENERATION {generation_n}==')
        print('selecting parents...')
        parents = select_parents(ants, n_parents, 15)

        print('generating children...')
        children = get_children(parents, probability_crossover, n_children)

        elites = select_elite(ants, fitness, n_elites)
        print('mutating children...')
        children, n_mutations = mutate(children, probability_mutation)
        print(f'{n_mutations} mutations occured!')

        assert(len(children + elites) == population_size)
        ants = children + elites
        print('running simulations...')
        for ant in ants:
            trial, fitness[ant], _ = ant_simulator(food_map, map_size, ant)

        ants = sorted(ants, key=lambda x: fitness[x], reverse=True)
        print(f'fitness statistics of final population after generation {generation_n}:')
        print('high:', fitness[ants[0]])
        print('75%:', fitness[ants[population_size//4]])
        print('50%:', fitness[ants[population_size//2]])
        print('25%:', fitness[ants[3*population_size//4]])
        print('low:', fitness[ants[population_size-1]])
        analytics['best_fitness_per_gen'].append(fitness[ants[0]])
        print()

    analytics['ants'] = ants

    print('==BEST ANT==')
    best = ants[0]

    #create heatmap for best ant
    trial, best_fitness, heatmap = ant_simulator(food_map, map_size, best, heatmap=True)
    analytics['heatmap'] = heatmap
    display_trials(trial, f"final.txt")
    display_trial_img(heatmap, map_size, f'final.png')
    print('heatmap saved to final.png')
    
    print('genes -- fitness')
    print(best, '--', best_fitness)
    
    
    return best, fitness, analytics

def initialize_population(num_population, greedy_proportion=0.1):
    ants = []
    n_greedy = round(greedy_proportion*num_population)
    for _ in range(n_greedy):
        # ant = ''
        # initialize ants greedily: state 0 moves forward, and 6 states redirect to state 0 if there is food directly ahead
        # more discussion in writeup
        ant = list('1x0' + 'ax0' * 5 + 'axx' * 4)
        for idx, c in enumerate(ant):
            # print(idx, c)
            if c == 'a':
                ant[idx] = str(random.randint(1, 4))
            if c == 'x':
                ant[idx] = str(random.randint(0, 9))
        ants.append(''.join(ant))
    for _ in range(num_population - n_greedy):
        # ant = ''
        # initialize some ants completely randomly
        ant = list('axx' * 10)
        for idx, c in enumerate(ant):
            # print(idx, c)
            if c == 'a':
                ant[idx] = str(random.randint(1, 4))
            if c == 'x':
                ant[idx] = str(random.randint(0, 9))
        ants.append(''.join(ant))
    return ants
    
    
def ant_simulator(food_map, map_size, ant_genes, heatmap=False):
    """
    parameters:
        food_map: a list of list of strings, representing the map of the environment with food
            "1": there is a food at the position
            "0": there is no food at the position
            (0, 0) position: the top left corner of the map
            (x, y) position: x is the row, and y is the column
        map_size: a list of int, the dimension of the map. It is in the format of [row, column]
        ant_genes: a string with length 30. It encodes the ant's genes, for more information, please refer to the handout.
    
    return:
        trial: a list of list of strings, representing the trials
            1: there is food at that position, and the spot was not visited by the ant
            0: there is no food at that position, and the spot was not visited by the ant
            empty: the spot has been visited by the ant
	fitness: fitness of ant
    heatmap: a list of list of 3-tuples, an RGB array for a heatmap of the trial
    
    It takes in the food_map and its dimension of the map and the ant's gene information, and return the trial in the map
    """
    color = None
    step_time = 200
    
    trial = []
    for i in food_map:
        line = []
        for j in i:
            line.append(j)
        trial.append(line)

    if heatmap:
        color = []
        for i in food_map:
            line = []
            for j in i:
                if j == '0':
                    line.append((0,0,0))
                else:
                    line.append((255, 40, 0))
            color.append(line)

    position_x, position_y = 0, 0
    orientation = [(1, 0), (0, -1), (-1, 0), (0, 1)] # face down, left, up, right
    fitness = 0
    state = 0
    orientation_state = 3
    gene_list = [ant_genes[i : i + 3] for i in range(0, len(ant_genes), 3)]
    
    for i in range(step_time):
        if heatmap:
            if trial[position_x][position_y] != "1":
                old = color[position_x][position_y]
                new = (old[0] + 40, old[1], old[2] + 40)
                blueness_inverse = 1 - i/400
                new = (round(new[0]*blueness_inverse), new[1], new[2])
            else:
                fitness += 1
                new = (0, 150, 0)
            color[position_x][position_y] = new
        else:
            if trial[position_x][position_y] == "1":
                fitness += 1

        trial[position_x][position_y] = " "
        
        sensor_x = (position_x + orientation[orientation_state][0]) % map_size[0]
        sensor_y = (position_y + orientation[orientation_state][1]) % map_size[1]
        sensor_result = trial[sensor_x][sensor_y]
        
        if sensor_result == "1":
            state = int(gene_list[state][2])
        else:
            state = int(gene_list[state][1])
        
        action = gene_list[state][0]
        
        if action == "1": # move forward
            position_x = (position_x + orientation[orientation_state][0]) % map_size[0]
            position_y = (position_y + orientation[orientation_state][1]) % map_size[1]
        elif action == "2": # turn right
            orientation_state = (orientation_state + 1) % 4
        elif action == "3": # turn left
            orientation_state = (orientation_state - 1) % 4
        elif action == "4": # do nothing
            pass
        else:
            raise Exception("invalid action number!")
    
    return trial, fitness, color
        

def get_map(file_name):
    """
    parameters:
        file_name: a string, the name of the file which stored the map. The first line of the map is the dimension (row, column), the rest is the map
            1: there is a food at the position
            0: there is no food at the position
    
    return:
        food_map: a list of list of strings, representing the map of the environment with food
            "1": there is a food at the position
            "0": there is no food at the position
            (0, 0) position: the top left corner of the map
            (x, y) position: x is the row, and y is the column
        map_size: a list of int, the dimension of the map. It is in the format of [row, column]
    
    It takes in the file_name of the map, and return the food_map and the dimension map_size
    """
    food_map = []
    map_file = open(file_name, "r")
    first_line = True
    map_size = []
    
    for line in map_file:
        line = line.strip()
        if first_line:
            first_line = False
            map_size = line.split()
            continue
        if line:
            food_map.append(line.split())
    
    map_file.close()
    return food_map, [int(i) for i in map_size]

def display_trials(trials, target_file):
    """
    parameters:
        trials: a list of list of strings, representing the trials
            1: there is food at that position, and the spot was not visited by the ant
            0: there is no food at that position, and the spot was not visited by the ant
            empty: the spot has been visited by the ant
        taret_file: a string, the name the target_file to be saved
    
    It takes in the trials, and target_file, and saved the trials in the target_file. You can open the target_file to take a look at the ant's trial.
    """
    trial_file = open(target_file, "w")
    for line in trials:
        trial_file.write(" ".join(line))
        trial_file.write("\n")
    trial_file.close()

def display_trial_img(colors, dims, filename):
    img = Image.new('RGB', dims[::-1])
    data = [pixel for row in colors for pixel in row]
    # print(data)
    img.putdata(data)
    img.save(filename)
    

if __name__ == "__main__":
    import sys
    trail = 'muir'
    population_size = 100
    max_generation = 1000
    probability_crossover = 0.3
    probability_mutation = 0.1
    parent_proportion = 0.6
    n_children = 3
    greedy_proportion = 0.1
    seed = 0

    if len(sys.argv) in (9, 10):
        trail = sys.argv[1]
        population_size = int(sys.argv[2])
        max_generation = int(sys.argv[3])
        probability_crossover = float(sys.argv[4])
        probability_mutation = float(sys.argv[5])
        parent_proportion = float(sys.argv[6])
        n_children = int(sys.argv[7])
        greedy_proportion = float(sys.argv[8])
        if len(sys.argv) == 10:
            seed = int(sys.argv[9])
    elif len(sys.argv) != 1:
        print(f'usage: {sys.argv[0]} [<trail> <n> <g> <c> <m> <pp> <nc> <gp> [<seed>]] ')
        exit(1)

    random.seed(seed)
    
    # run GA
    food_map, map_size = get_map(f"{trail}.txt")
    best, fitness, analytics = genetic_algorithm(population_size, food_map, map_size, max_generation, probability_crossover, 
        probability_mutation, parent_proportion, n_children, greedy_proportion)

    # print report
    from pathlib import Path
    ants = analytics['ants']
    file_prefix = f'{trail}/n{population_size}-g{max_generation}-c{probability_crossover}-m{probability_mutation}-p{parent_proportion}-nc{n_children}-gp{greedy_proportion}-s{seed}'
    Path(f'reports/{file_prefix}').mkdir(parents=True, exist_ok=True)
    with open(f'reports/{file_prefix}/report.txt', 'w') as f:
        print('==params==', file=f)
        print('population_size:', population_size, file=f)
        print('max_generation:', max_generation, file=f)
        print('probaility_crossover:', probability_crossover, file=f)
        print('probability_mutation:', probability_mutation, file=f)
        print('parent_proportion:', parent_proportion, file=f)
        print('n_children:', n_children, file=f)
        print('greedy_proportion:', greedy_proportion, file=f)
        print(file=f)
        print('==results==', file=f)
        print('high:', fitness[ants[0]], file=f)
        print('75%:', fitness[ants[population_size//4]], file=f)
        print('50%:', fitness[ants[population_size//2]], file=f)
        print('25%:', fitness[ants[3*population_size//4]], file=f)
        print('low:', fitness[ants[population_size-1]], file=f)
        print(file=f)
        print('==best ant genes==', file=f)
        print(ants[0], file=f)
    
    display_trial_img(analytics['heatmap'], map_size, f'reports/{file_prefix}/heatmap.png')

    # plot fitness over generation
    gens = list(range(1, max_generation + 1))
    fitness_ = analytics['best_fitness_per_gen']
    fig, ax = plt.subplots()
    ax.plot(gens, fitness_)

    ax.set(xlabel='generation number', ylabel='fitness of best ant',
        title=f'fitness over generation {file_prefix}')
    ax.grid()

    fig.savefig(f"reports/{file_prefix}/graph.png")
    # plt.show()
    print('report saved under', file_prefix)

    #comment this out to generate the final generation plots
    exit()


    #run santafe trials
    fitness_2 = []
    food_map_santafe, map_size_santafe = get_map("santafe.txt")
    for i, ant in enumerate(ants):
        trial, f, heatmap = ant_simulator(food_map_santafe, map_size_santafe, ant, heatmap=True)
        display_trials(trial, f"santafe/santafe-{i}.txt")
        display_trial_img(heatmap, map_size_santafe, f'santafe/santafe-{i}.png')
        fitness_2.append(f)
    


    #plot fitness of final gen for both trials
    fitness_ = [fitness[ant] for ant in ants]

    fig, ax = plt.subplots()
    ax.bar(list(range(len(ants))), fitness_)
    plt.xticks([0,5,10,15,20,25,30,35,40,45,50])
    plt.yticks(range(0, 90, 5))
    ax.set(xlabel='ant number', ylabel='fitness',
        title='Fitness of final ant population on Muir trial')
    ax.grid(which='both')

    fig.savefig("fitness_final_gen.png")
    plt.show()
    

    fig, ax = plt.subplots()
    ax.bar(list(range(len(ants))), fitness_2)
    plt.xticks([0,5,10,15,20,25,30,35,40,45,50])
    plt.yticks(range(0, 90, 5))
    ax.set(xlabel='ant number', ylabel='fitness',
        title='Fitness of final ant population on Santa Fe trial')
    ax.grid(which='both')

    fig.savefig("fitness_final_gen_santafe.png")
    plt.show()