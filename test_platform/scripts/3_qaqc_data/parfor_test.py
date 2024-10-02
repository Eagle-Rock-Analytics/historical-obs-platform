# import the parfor function; note
# that this will automatically initialize MPI
# also import pprint for parallel-friendly printing
from simplempi.parfor import parfor, pprint
import time

# define a list to loop over
my_list = list(range(100)) 

# define a function that does something with each item in my_list
def func(i):
    if i%2==0:
        time.sleep(2)
        return i**2
    else:
        time.sleep(2)
        return i**2

# loop in parallel over my_list
for i in parfor(my_list):
    result = func(i)
    pprint(f"{i}**2 = {result}", flush=True)
