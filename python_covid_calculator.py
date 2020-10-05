import math
import random
import matplotlib.pyplot as plt

parameters = {"breathing_rate":0.66,"duration_of_event":0.8,"fraction_people_with_masks":1,"inhalation_mask_efficiency":0.3,"width":9,"length":6,"height":3,"ventilations_per_hour":1,"decay_rate_of_virus":0.62,"deposition_to_surfaces":0.3,"additional_control_measures":0,"exhalation_mask_efficiency":0.5,"quanta_exhalation_rate":135}

def initialize_params(params):
    total_first_order_loss_rate = params["ventilations_per_hour"] + params["decay_rate_of_virus"] + params["deposition_to_surfaces"] + params["additional_control_measures"]
    volume = params["width"] * params["length"] * params["height"]
    return total_first_order_loss_rate,volume,params["quanta_exhalation_rate"],params["exhalation_mask_efficiency"],params["fraction_people_with_masks"],params["duration_of_event"],params["breathing_rate"],params["inhalation_mask_efficiency"]

def calculate_infection_rate (susceptible_people,infective_people):
    global total_first_order_loss_rate,volume,quanta_exhalation_rate,exhalation_mask_efficiency,fraction_people_with_masks,duration_of_event,breathing_rate,inhalation_mask_efficiency
    net_emission_rate = quanta_exhalation_rate * (1 - exhalation_mask_efficiency * fraction_people_with_masks) * infective_people
    avg_quanta_concentration = (net_emission_rate / total_first_order_loss_rate / volume) * (1 - (1 /total_first_order_loss_rate / duration_of_event) * (1 - math.exp(-total_first_order_loss_rate * duration_of_event)))
    quanta_inhaled_per_person = avg_quanta_concentration * breathing_rate * duration_of_event * (1 - inhalation_mask_efficiency * fraction_people_with_masks)
    infection_likelihood = 1 - math.exp(-quanta_inhaled_per_person)
    return infection_likelihood * susceptible_people

def conditional_likelihood_evolution():
    days = 1
    infected = [1]
    infection_rate = [0]
    new_infected = [0]
    susceptible_people = 15
    remaining_likelihood = 0
    infective_people = 1
    while (susceptible_people > 0) :
        if (days > 14):
            susceptible_people -= infected[days - 15] - new_infected[days - 15]
        infection_rate_for_the_day = calculate_infection_rate(susceptible_people,infective_people)
        remaining_likelihood,new_daily_infected = math.modf(infection_rate_for_the_day + remaining_likelihood)
        infected.append(infected[days - 1] + new_daily_infected)
        infection_rate.append(infection_rate_for_the_day)
        new_infected.append(new_daily_infected)
        susceptible_people = susceptible_people - new_daily_infected
        infective_people = infective_people + new_daily_infected
        days += 1
        print("Day {} - Infected : {} - Susceptible people : {} - Infection rate : {}".format(days - 1,infected[days - 1],susceptible_people,infection_rate_for_the_day))

def calculate_newly_infected(susceptible_people,infective_people):
    global total_first_order_loss_rate,volume,quanta_exhalation_rate,exhalation_mask_efficiency,fraction_people_with_masks,duration_of_event,breathing_rate,inhalation_mask_efficiency
    net_emission_rate = quanta_exhalation_rate * (1 - exhalation_mask_efficiency * fraction_people_with_masks) * infective_people
    avg_quanta_concentration = (net_emission_rate / total_first_order_loss_rate / volume) * (1 - (1 /total_first_order_loss_rate / duration_of_event) * (1 - math.exp(-total_first_order_loss_rate * duration_of_event)))
    quanta_inhaled_per_person = avg_quanta_concentration * breathing_rate * duration_of_event * (1 - inhalation_mask_efficiency * fraction_people_with_masks)
    infection_likelihood = 1 - math.exp(-quanta_inhaled_per_person)
    newly_infected = 0
    for i in range(susceptible_people):
        likelihood = random.random()
        if likelihood < infection_likelihood:
            newly_infected += 1
    return newly_infected,infection_likelihood

def simulate(iterations):
    days_simulated = []
    for i in range(iterations):
        days_simulated.append(likelihood_evolution())
    arithmetic_mean = sum(days_simulated)*1.0/len(days_simulated)
    print("DAYS : "+str(arithmetic_mean))

def likelihood_evolution():
    days = 1
    infected = [0]
    infection_rate = [0]
    new_infected = [0]
    susceptible_people = 15
    infective_people = 1
    while (susceptible_people > 0) :
        if (days > 14):
            susceptible_people -= infected[days - 15] - new_infected[days - 15]
        new_daily_infected,infection_likelihood = calculate_newly_infected(susceptible_people,infective_people)
        infected.append(infected[days - 1] + new_daily_infected)
        new_infected.append(new_daily_infected)
        susceptible_people = susceptible_people - new_daily_infected
        infective_people = infective_people + new_daily_infected
        days += 1
        #print("Day {} - Infected : {} - Susceptible people : {} - Likelihood : {}".format(days - 1,infected[days - 1],susceptible_people,infection_likelihood))
    return days


total_first_order_loss_rate,volume,quanta_exhalation_rate,exhalation_mask_efficiency,fraction_people_with_masks,duration_of_event,breathing_rate,inhalation_mask_efficiency = initialize_params(parameters)
conditional_likelihood_evolution()
print("===============================================================================================")
#likelihood_evolution()
simulate(1000000)
#plt.plot(range(days),infected,'b',label='Infected')
#plt.plot(range(days),infection_rate,'r',label='Infection rate')
#plt.legend(loc='upper left')
#plt.xlabel('Days')
#plt.ylabel('Rate')
#plt.grid(True)
#plt.show()
