from scipy import stats
from scipy.stats import t as t_dist
from scipy.stats import chi2

from abtesting_test import *

##########################
# NOTE TO TA'S: I used Repl.it for this assignment, so if there happens to be some discrepency of 
# formatting/importing/testing on this file etc, it is likely due to that. I've copied my code from my Repl.it repl 
# and pasted it onto this python file, and decided to include everything so you may see my complete process and tests,
# which is why I have print() lines (which is how I ran tests on Repl.it) for my tests at the bottom! 
##########################

def slice_2D(list_2D, start_row, end_row, start_col, end_col):
    '''
    Splices a the 2D list via start_row:end_row and start_col:end_col
    :param list: list of list of numbers
    :param nums: start_row, end_row, start_col, end_col
    :return: the spliced 2D list (ending indices are exclsive)
    '''
    to_append = []
    for l in range(start_row, end_row):
        to_append.append(list_2D[l][start_col:end_col])

    return to_append

def get_avg(nums):
    '''
    Helper function for calculating the average of a sample.
    :param nums: list of numbers
    :return: average of list
    '''
    sum = 0
    for i in nums:
      sum += i
    avg = sum / len(nums)
    return avg


def get_stdev(nums):
    '''
    Helper function for calculating the standard deviation of a sample.
    :param nums: list of numbers
    :return: standard deviation of list
    '''
    dev_sum = 0
    avg = get_avg(nums)
    for i in nums:
      dev_sum += (i - avg)**2
    
    return (dev_sum/(len(nums)-1))**(0.5)


def get_standard_error(a, b):
    '''
    Helper function for calculating the standard error, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: standard error of a and b (see studio 6 guide for this equation!)
    '''
    std_a = get_stdev(a)
    std_b = get_stdev(b)
    n_a = len(a)
    n_b = len(b)
    return ((((std_a)**2) / n_a) + (((std_b)**2) / n_b))**(1/2)


def get_2_sample_df(a, b):
    '''
    Calculates the combined degrees of freedom between two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: integer representing the degrees of freedom between a and b (see studio 6 guide for this equation!)
    HINT: you can use Math.round() to help you round!
    '''
    se = get_standard_error(a,b)
    std_a = get_stdev(a)
    std_b = get_stdev(b)
    n_a = len(a)
    n_b = len(b)
    div_a = (((std_a**2)/n_a)**2)/(n_a - 1)
    div_b = (((std_b**2)/n_b)**2)/(n_b - 1)
    freedom = (se**4)/(div_a + div_b)
    return round(freedom)
    

def get_t_score(a, b):
    '''
    Calculates the t-score, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: number representing the t-score given lists a and b (see studio 6 guide for this equation!)
    '''
    se = get_standard_error(a,b)
    t = (get_avg(a) - get_avg(b))/se
    if(t > 0):
      t = t * -1
    return t

def perform_2_sample_t_test(a, b):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates a p-value by performing a 2-sample t-test, given two lists of numbers.
    :param a: list of numbers
    :param b: list of numbers
    :return: calculated p-value
    HINT: the t_dist.cdf() function might come in handy!
    '''
    t_score_n = get_t_score(a,b)
    df = get_2_sample_df(a,b)

    return t_dist.cdf(t_score_n, df)
    
# [OPTIONAL] Some helper functions that might be helpful in get_expected_grid().
def row_sum(observed_grid, ele_row):
  r_sum = 0
  for i in range(len(observed_grid[0])):
    r_sum += observed_grid[ele_row][i]
  return r_sum

def col_sum(observed_grid, ele_col):
  c_sum = 0
  for i in range(0, len(observed_grid)):
    c_sum += observed_grid[i][ele_col]
  return c_sum

def total_sum(observed_grid):
  t_sum = 0
  for i in range(0, len(observed_grid)):
    for j in range(0, len(observed_grid[0])):
      t_sum += observed_grid[i][j]

  return t_sum

def calculate_expected(row_sum, col_sum, tot_sum):
  return (row_sum * col_sum) / tot_sum

def get_expected_grid(observed_grid):
    '''
    Calculates the expected counts, given the observed counts.
    ** DO NOT modify the parameter, observed_grid. **
    :param observed_grid: 2D list of observed counts
    :return: 2D list of expected counts
    HINT: To clean up this calculation, consider filling in the optional helper functions below!
    '''
    
    filled_grid = [[0] * len(observed_grid[0]) for i in range(len(observed_grid))]

    for i in range(len(observed_grid)):
      for j in range(len(observed_grid[0])):
        filled_grid[i][j] = calculate_expected(row_sum(observed_grid, i), col_sum(observed_grid, j), total_sum(observed_grid))

    return filled_grid

def df_chi2(observed_grid):
    '''
    Calculates the degrees of freedom of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: degrees of freedom of expected counts (see studio 6 guide for this equation!)
    '''
    rows = len(observed_grid)
    columns = len(observed_grid[0])
    return (rows-1)*(columns-1)

def chi2_value(observed_grid):
    '''
    Calculates the chi^2 value of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: associated chi^2 value of expected counts (see studio 6 guide for this equation!)
    '''
    expected = get_expected_grid(observed_grid)
    chi = 0

    for i in range(len(observed_grid)):
      for j in range(len(observed_grid[0])):
        chi += (((observed_grid[i][j] - expected[i][j])**2) / expected[i][j])
    
    return chi

def perform_chi2_homogeneity_test(observed_grid):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates the p-value by performing a chi^2 test, given a list of observed counts
    :param observed_grid: 2D list of observed counts
    :return: calculated p-value
    HINT: the chi2.cdf() function might come in handy!
    '''
    return round(1-chi2.cdf(chi2_value(observed_grid), df_chi2(observed_grid)), 6)

def data_to_num_list(s):
  '''
    Takes a copy and pasted row/col from a spreadsheet and produces a usable list of nums. 
    This will be useful when you need to run your tests on your cleaned log data!
    :param str: string holding data
    :return: the spliced list of numbers
    '''
  return list(map(float, s.split()))


#######
# EXAMPLE TESTS PROVIDED BY STENCIL CODE
#######

# t_test 1:
a_t1_list = data_to_num_list(a1)
b_t1_list = data_to_num_list(b1)
print(((get_t_score(a_t1_list, b_t1_list)))) # this should be -129.500
print(perform_2_sample_t_test(a_t1_list, b_t1_list)) # this should be 0.0000
# why do you think this is? Take a peek at a1 and b1 in abtesting_test.py :)

# t_test 2:
a_t2_list = data_to_num_list(a2)
b_t2_list = data_to_num_list(b2)
print(get_t_score(a_t2_list, b_t2_list)) # this should be -1.48834
print(perform_2_sample_t_test(a_t2_list, b_t2_list)) # this should be .082379

# t_test 3:
a_t3_list = data_to_num_list(a3)
b_t3_list = data_to_num_list(b3)
print(get_t_score(a_t3_list, b_t3_list)) # this should be -2.88969
print(perform_2_sample_t_test(a_t3_list, b_t3_list)) # this should be .005091

# chi2_test 1:
a_c1_list = data_to_num_list(a_count_1) 
b_c1_list = data_to_num_list(b_count_1)
c1_observed_grid = [a_c1_list, b_c1_list]
print(chi2_value(c1_observed_grid)) # this should be 4.103536
print(perform_chi2_homogeneity_test(c1_observed_grid)) # this should be .0427939

# chi2_test 2:
a_c2_list = data_to_num_list(a_count_2) 
b_c2_list = data_to_num_list(b_count_2)
c2_observed_grid = [a_c2_list, b_c2_list]
print(chi2_value(c2_observed_grid)) # this should be 33.86444
print(perform_chi2_homogeneity_test(c2_observed_grid)) # this should be 0.0000
# Again, why do you think this is? Take a peek at a_count_2 and b_count_2 in abtesting_test.py :)

# chi2_test 3:
a_c3_list = data_to_num_list(a_count_3) 
b_c3_list = data_to_num_list(b_count_3)
c3_observed_grid = [a_c3_list, b_c3_list]
print(chi2_value(c3_observed_grid)) # this should be .3119402
print(perform_chi2_homogeneity_test(c3_observed_grid)) # this should be .57649202


####################
# MY COLLECTED DATA
####################

a_completed = """
13261
5275
14974
7260
0
"""
b_completed = """
4207
9651
23146
14963
10940
8428
18027
8012
"""
a_return = """
1	4
"""
b_return = """
4	4
"""

######
# RUNNING MY OWN DATA THROUGH THE TESTS
######

a_my_list = data_to_num_list(a_completed) 
b_my_list = data_to_num_list(b_completed)
print(get_t_score(a_my_list, b_my_list))
print(perform_2_sample_t_test(a_my_list, b_my_list))


a_my_return = data_to_num_list(a_return)
b_my_return = data_to_num_list(b_return)
my_observed_grid = [a_my_return, b_my_return]
print(chi2_value(my_observed_grid)) 
print(perform_chi2_homogeneity_test(my_observed_grid))