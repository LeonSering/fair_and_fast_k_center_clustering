from ff_k_center import FFKCenter
import numpy as np

model = FFKCenter(4,1,[(1,13)]) # parameters: k, privacy_bound = 1, rep_intervals = []
pos = np.array([[1,2],[0,0],[6,7],[0,-1],[-2,-3],[0,-1.2],[1,3],[6,1]])
colors = np.array([0,0,0,1,1,0,0,1])


model.fit(pos, colors, verbose = 2, thread_count = 6, phase_2_rerun = False, phase_5_gonzalez = True) # positions of points; followed by colors; Optional: verbose (0: silent, 1: brief, 2: verbose; default: 1); thread_count (default: #cores); phase_2_rerun (default: True); phase_5_gonzalez (default: False);

print("\nr^*: ", model.radius) # the final radius of the assignment
print("C^*: ", model.centers) # Chosen centers by the point-index.
print("phi^*: ", model.assignment) # for each point the center it is assigned to.
print("running time: ", model.running_time) # running time in sec
