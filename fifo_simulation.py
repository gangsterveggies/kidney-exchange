import numpy as np
import random as rand

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

  @staticmethod  
  def sample_type():
    return utils.weighted_choice(blood_type.perc_list)

class dp_pair:
  def __init__(self, _donor, _patient):
    self.donor = _donor
    self.patient = _patient

  def is_compatible(self, other_pair):
    return self.donor.is_patient_compatible(other_pair.patient) and self.patient.is_donor_compatible(other_pair.donor)

if __name__ == "__main__":
  print blood_type.sample_type()
  print utils.sample_exponential()
