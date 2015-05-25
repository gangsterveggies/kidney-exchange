import numpy as np
import random as rand
import scipy.stats as scist
import networkx as nx

_O = 0
_A = 1
_B = 2
_AB = 3

class utils:
  @staticmethod
  def weighted_choice(choices):
    r = np.random.uniform()
    upto = 0
    for c, w in choices:
      if upto + w > r:
        return c
      upto += w

  @staticmethod
  def sample_exponential(avg=2):
    return np.random.exponential(1.0/avg)

  @staticmethod
  def mean_confidence_interval(data, alpha):
    n = len(data)
    m = float(sum(data)) / n
    var = sum([(x - m)**2 for x in data]) / float(n - 1)

    tfact = scist.t._ppf(1 - alpha / 2., n - 1)
    h = tfact * np.sqrt(var / n)

    return (m - h, m + h)

  @staticmethod
  def padded_print(initial_text, array, pad):
    print initial_text,
    for i in array:
      print str(i).center(pad, " "),
    print ""

  @staticmethod
  def make_match(pair_list, stats):
    G = nx.Graph()
    G.add_nodes_from(range(len(pair_list)))

    for i in range(len(pair_list)):
      for j in range(i + 1, len(pair_list)):
        if pair_list[i].is_compatible(pair_list[j]):
          G.add_edge(i, j)

    matches = nx.max_weight_matching(G)
    matches_pairs = matches.items()
    matched = 0

    for i in range(len(matches_pairs)):
      if matches_pairs[i][0] <= matches_pairs[i][1]:
        if matched < len(matches_pairs) / 2 - 1:
          stats.exchange(pair_list[matches_pairs[i][0]], pair_list[matches_pairs[i][1]])
        else:
          stats.exchange(pair_list[matches_pairs[i][0]], pair_list[matches_pairs[i][1]], change_size=True)
        matched += 1

    return [pair_list[i] for i in range(len(pair_list)) if i not in matches]

class blood_type:
  comp_matrix = [[True, True, True, True], [False, True, False, True], [False, False, True, True], [False, False, False, True]]
  perc_list = [(_O, 0.423), (_A, 0.464), (_B, 0.078), (_AB, 0.035)]

  def __init__(self, _btype):
    if (_btype >= 4 or _btype < 0):
      raise IndexError
    self.btype = _btype

  def is_patient_compatible(self, patient):
    return blood_type.comp_matrix[self.btype][patient.btype]

  def is_donor_compatible(self, donor):
    return blood_type.comp_matrix[donor.btype][self.btype]

  def __str__(self):
    return ['O', 'A', 'B', 'AB'][self.btype]

  @staticmethod  
  def sample_type():
    return blood_type(utils.weighted_choice(blood_type.perc_list))

class dp_pair:
  def __init__(self, _donor, _patient):
    self.donor = _donor
    self.patient = _patient

  def is_compatible(self, other_pair):
    return self.donor.is_patient_compatible(other_pair.patient) and self.patient.is_donor_compatible(other_pair.donor)

  def __str__(self):
    return "(D: %s, P: %s)" % (str(self.donor), str(self.patient))

  @staticmethod
  def sample_dp_pair():
    npair = dp_pair(blood_type.sample_type(), blood_type.sample_type())
    next_time = utils.sample_exponential()
    return (npair, next_time)

class stats_recorder:
  def __init__(self):
    self.avg_size = 0
    self.current_time = 0
    self.current_count = 0
    self.total_time = 0
    self.donor_count = [0, 0, 0, 0]
    self.patient_count = [0, 0, 0, 0]
    self.exchange_matrix = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    self.total_exchanges = 0

  def enter(self, pair, next_time):
    self.donor_count[pair.donor.btype] += 1
    self.patient_count[pair.patient.btype] += 1
    self.current_time += next_time
    self.total_time += next_time
    self.current_count += 1

  def exchange(self, pair1, pair2, self_exchange=False, change_size=False):
    if change_size:
      self.avg_size += self.current_count * self.current_time
    self.exchange_matrix[pair1.donor.btype][pair2.patient.btype] += 1
    if not(self_exchange):
      self.exchange_matrix[pair2.donor.btype][pair1.patient.btype] += 1

    if self_exchange:
      self.current_count -= 1
    else:
      self.current_count -= 2

    if change_size:
      self.current_time = 0

  def finalize(self):
    self.avg_size += self.current_count * self.current_time    
    self.avg_size /= self.total_time
    self.current_time = 0

def simulate(total_time=156):
  time = 0
  next_match_time = 4
  waiting_queue = []
  stats = stats_recorder()

  while time < total_time:
    if time >= next_match_time:
      waiting_queue = utils.make_match(waiting_queue, stats)
      next_match_time += 4

      for i in range(len(waiting_queue) - 1, -1, -1):
        if waiting_queue[i].is_compatible(waiting_queue[i]):
          stats.exchange(waiting_queue[i], waiting_queue[i], self_exchange=True, change_size=True)
          del waiting_queue[i]

    npair, next_time = dp_pair.sample_dp_pair()

    stats.enter(npair, next_time)
    waiting_queue.append(npair)
    time += next_time

  if total_time >= next_match_time:
    waiting_queue = utils.make_match(waiting_queue, stats)
    next_match_time += 4

    for i in range(len(waiting_queue) - 1, -1, -1):
      if waiting_queue[i].is_compatible(waiting_queue[i]):
        stats.exchange(waiting_queue[i], waiting_queue[i], self_exchange=True, change_size=True)
        del waiting_queue[i]

  stats.finalize()

  return (stats.avg_size, stats.current_count, stats.donor_count, stats.patient_count, stats.exchange_matrix)

if __name__ == "__main__":
  max_iter = 100
  avg_size_list = []
  final_list = []
  donor_list = []
  patient_list = []
  exchange_list = []

  for i in range(max_iter):
    avg, final, donor, patient, exchange = simulate()

    avg_size_list.append(avg)
    final_list.append(final)
    donor_list.append(np.array(donor))
    patient_list.append(np.array(patient))
    exchange_list.append(np.array(exchange))

  print "        Homework 3 - Problem 3 / alternative"
  print "     Joao Ramos, Pedro Paredes"
  print ""

  print "Showing results after %d simulations" % max_iter
  print "Simulated over 3 years"
  print ""

  print "Average number of pairs in the waiting queue per week: %d" % round(sum(avg_size_list) / float(max_iter))
  print "95%% Confidence interval for the average number of pairs in the waiting list per week: [%0.2f, %0.2f]" % utils.mean_confidence_interval(avg_size_list, 0.05)
  print ""

  print "Average number of pairs in the waiting queue after the 3 years: %d" % round(sum(final_list) / float(max_iter))
  print "95%% Confidence interval for the average number of pairs in the waiting list after the 3 years: [%0.2f, %0.2f]" % utils.mean_confidence_interval(final_list, 0.05)
  print ""

  print "Donor by blood type count average:"
  utils.padded_print("Types:", ['O', 'A', 'B', 'AB'], 7)
  utils.padded_print("      ", np.round(sum(donor_list) / float(max_iter)), 7)
  print ""

  print "Patient by blood type count average:"
  utils.padded_print("Types:", ['O', 'A', 'B', 'AB'], 7)
  utils.padded_print("      ", np.round(sum(patient_list) / float(max_iter)), 7)
  print ""

  print "Average exchanges by blood type in donor x patient matrix:"
  avg_ex = sum(exchange_list) / float(max_iter)
  print "                    Patient"
  utils.padded_print("Donor:", ['O', 'A', 'B', 'AB'], 7)
  utils.padded_print("  O  |", np.round(avg_ex[0], 2), 7)
  utils.padded_print("  A  |", np.round(avg_ex[1], 2), 7)
  utils.padded_print("  B  |", np.round(avg_ex[2], 2), 7)
  utils.padded_print("  AB |", np.round(avg_ex[3], 2), 7)
  print ""

  avg_total = sum(sum(sum(exchange_list) / float(max_iter)))
  print "Total average number of exchanges: %0.2f" % avg_total
  print ""

  print "Normalized average exchanges by blood type in donor x patient matrix:"
  avg_ex /= avg_total
  print "                    Patient"
  utils.padded_print("Donor:", ['O', 'A', 'B', 'AB'], 7)
  utils.padded_print("  O  |", np.round(avg_ex[0], 2), 7)
  utils.padded_print("  A  |", np.round(avg_ex[1], 2), 7)
  utils.padded_print("  B  |", np.round(avg_ex[2], 2), 7)
  utils.padded_print("  AB |", np.round(avg_ex[3], 2), 7)
  print ""
