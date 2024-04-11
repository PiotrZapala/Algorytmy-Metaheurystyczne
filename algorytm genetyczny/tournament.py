def tournament_selection(self):
        competition1 = []
        competition2 = []
        individuals1 = []
        individuals2 = []
        total_fitness1 = 0
        total_fitness2 = 0
        size_of_population = len(self.set_of_individuals)
        values = list(range(size_of_population))
        for _ in range(int(size_of_population/2)):
            for _ in range(5):
                rand_father = choice(values)
                total_fitness1 = total_fitness1 + self.set_of_individuals[rand_father].fitness
                competition1.append(self.set_of_individuals[rand_father].fitness)
                individuals1.append(self.set_of_individuals[rand_father])
                rand_mother = choice(values)
                total_fitness2 = total_fitness2 + self.set_of_individuals[rand_mother].fitness
                competition2.append(self.set_of_individuals[rand_mother].fitness)
                individuals2.append(self.set_of_individuals[rand_mother])
            for n in range(len(competition1)):
                self.probability.append(competition1[n]/total_fitness1)
                self.field = self.field + self.probability[n]
                self.roulette.append(self.field)           
            r1 = uniform(0,1)
            for m in range(len(self.roulette)):
                if r1 <= self.roulette[m]:
                    self.selected_fathers.append(individuals1[m])
                    break
                else:
                    continue  
            self.probability = []
            self.field = 0
            self.roulette = []
            for o in range(len(competition2)):
                self.probability.append(competition2[o]/total_fitness2)
                self.field = self.field + self.probability[o]
                self.roulette.append(self.field)
            r2 = uniform(0,1)
            for p in range(len(self.roulette)):
                if r2 <= self.roulette[p]:
                    self.selected_mothers.append(individuals2[p])
                    break
                else:
                    continue
            self.probability = []
            self.field = 0
            self.roulette = []

        return self.selected_fathers, self.selected_mothers