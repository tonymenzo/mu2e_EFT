from num_input import NumericalInput
import wilson_coefficients as wc


# Define the default input parameters
default_input = NumericalInput().input_parameters

# Wrapper class
class WC_3f(wc.WC_3_flavor):
    def __init__(self, coeff_dict, user_input_dict=None):
        """ 
        'wrapper' class providing input for three-flavor Wilson coefficients 
        """
        if user_input_dict is None:
            self.ip = default_input
        else:
            self.ip = NumericalInput(user_input_dict).input_parameters
        wc.WC_3_flavor.__init__(self, coeff_dict, self.ip)