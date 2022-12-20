from pandas import read_csv
from gurobipy import Model, quicksum, GRB
import collections
import sys


# Authors: Mr Liam Timms, Dr Michael Forbes

# Change log: Version 8
# Modified: 29 Sep 2021
# Prev. version - CORPS_Inventory_V7_rearrange.py
# Allow program to take command line arguments to use batch file for cluster

# Function to take data filename as input, output file data as dictionary
def read_file(absolute, folder, name):
    # abs (string): absolute location of folder
    # folder (string): name of folder with data set
    # name (string): file name of data set without extension , e.g. 'A_sd'

    # Create list of variables in the file
    index = []
    for v in name.partition('_')[2]:
        index.append(v)

    # Read in csv file as dataframe
    df = read_csv(absolute + "\\" + folder + "\\" + name + '.csv', index_col=index)

    # Strip whitespace from column name
    df.columns = df.columns.str.replace(' ', '')

    # Return dictionary
    return df[name].to_dict()


# Function to convert duration (mins) to no. of time periods
def dur2period(num, period):
    # Inputs: num = number to round
    #         period = length of time period
    # Output:
    #         convert num to number of time periods (float)

    # Round num to nearest multiple of period
    nPeriods = round(num / period)

    return nPeriods


# Number of patients
# sizes = [35, 40, 45, 70, 80, 90, 105, 120, 135]size = sys.argv[1]
#
# trial = sys.argv[2]
#
# t_period = int(sys.argv[3])

# Trial number
# trials = range(1, 6)

# Length of time periods (for rounding durations)
# t_period = float(input('Enter time period length (mins): '))

# Absolute location of data file folder CORPS-Data\Alpha
# file_loc = r'C:\Users\admin\Documents\University\OR Project\CORPS-Data\Alpha'

# Program version no. (for file names)
version = '8'


# Optimise over data set 'size_i_trial_j', data located in absolute path file_loc
def surgery_optimise(size_no, trial_no, t_period, file_loc):
    """
    :param size_no: integer, number of patients (35, 40, 45, 70, 80, 90, 105, 120, or 135)
    :param trial_no: integer, trial number (1, 2, 3, 4, or 5)
    :param t_period: float, length of time periods
    :param file_loc: string, absolute path to data folder "CORPS-Data/Alpha", not including "CORPS-Data/Alpha"
    :return: optimise MIP problem, save log and solution to .txt file
    """
    file_loc += '\CORPS-Data\Alpha'

    # Name of subfolder with data set
    folder_name = 'size_' + str(size_no) + '_trial_' + str(trial_no)

    # Erase previous copy of output file
    with open('optimised_' + folder_name + '.txt', 'w+') as f:
        f.truncate(0)

    # Read in data files

    A_sd = read_file(file_loc, folder_name, 'A_sd')
    B_hdr = read_file(file_loc, folder_name, 'B_hdr')
    C_hdr = read_file(file_loc, folder_name, 'C_hdr')  # Note: C_hdr = 3/16 * K_hdr
    E_ps = read_file(file_loc, folder_name, 'E_ps')
    F_p = read_file(file_loc, folder_name, 'F_p')
    G_p = read_file(file_loc, folder_name, 'G_p')
    K_hdr = read_file(file_loc, folder_name, 'K_hdr')
    L_shd = read_file(file_loc, folder_name, 'L_shd')
    Omega_sp = read_file(file_loc, folder_name, 'Omega_sp')
    Q_phr = read_file(file_loc, folder_name, 'Q_phr')
    T_ps = read_file(file_loc, folder_name, 'T_ps')
    Theta_p = read_file(file_loc, folder_name, 'Theta_p')
    U_p = read_file(file_loc, folder_name, 'U_p')
    V_hdr = read_file(file_loc, folder_name, 'V_hdr')

    # Round durations E_ps and T_ps (same indices)

    for dur in E_ps:
        E_ps[dur] = dur2period(E_ps[dur], t_period)
        T_ps[dur] = dur2period(T_ps[dur], t_period)

    # Round prep time F_p and cleaning time G_p
    for p_index in F_p:
        F_p[p_index] = dur2period(F_p[p_index], t_period)
        G_p[p_index] = dur2period(G_p[p_index], t_period)

    Duration = [T_ps[p, s] for p, s in T_ps]

    minDuration = min(Duration)


    # print('Duration:', Duration)

    # --------------------------- SETS ----------------------------------------
    # Using indices from data sets
    # NOTE - Begin indexing at 1, e.g. days start from 1

    # Calculate no. of time periods
    DaySlots = list(set([B_hdr[k] for k in B_hdr]))  # Set up list of regular minutes for each hdr
    OvertimeSlots = list(set([V_hdr[k] for k in V_hdr]))  # Set up list of overtime minutes for each hdr
    if len(DaySlots) >= 1.1 or len(OvertimeSlots) >= 1.1:
        raise ValueError('Different available minutes for hdr')

    # Set B - regular time slots
    DaySlots = int(DaySlots[0] / t_period)  # Round regular into time slots
    # Set V - overtime slots
    OvertimeSlots = int(OvertimeSlots[0] / t_period)  # Round overtime into time slots
    # 8 hours regular time, 2 hours overtime
    nTime = DaySlots + OvertimeSlots  # number of time periods (t_period mins long)
    Times = range(1, nTime + 1)  # 8 hour day split into this many time periods

    # Patients
    P = range(1, size_no + 1)

    # Surgeons
    S = range(1, max([i[0] for i in A_sd.keys()]) + 1)

    # Hospitals
    H = range(1, max([i[0] for i in B_hdr.keys()]) + 1)

    # ORs in hospital h in H
    R_h = {h: list(set([i[2] for i in B_hdr.keys() if i[0] == h])) for h in H}

    # Days
    D = range(1, max([i[1] for i in A_sd.keys()]) + 1)

    # Patients that can be operated on in OR r of hospital h with theta_p >= d
    # Not necessary since ORs are homogenous in inventory model (except for cost)

    # Patients belonging to surgeon s and theta_p >= d, and surgeon s can operate on patient p
    # All surgeons can do all patients - implication of E_ps and T_ps
    Lambda_sd = {(s, d): [p for p in P if Theta_p[p] >= d and Omega_sp[s, p] == 1] for s in S for d in D}

    m = Model('OR Inventory Scheduling')
    m.setParam('Seed', 42)  # Set the Gurobi seed
    m.setParam('MIPGap', 0.0)  # Solve to optimality
    # m.setParam('BranchDir',1)
    # m.setParam('MIPFocus',2)
    m.setParam('TimeLimit', 21600)
    m.setParam('Threads', 4)
    # Save output log to same file as solution

    # Save output log to separate file to solution
    m.setParam('LogFile', 'V' + version + '_' + str(t_period) + 'min_' + folder_name + 'log.txt')

    # --------------------------- VARIABLES ----------------------------------------
    # Operation p is done by surgeon s in hospital h on day d starting at time t
    X = {(p, s, h, d, t): m.addVar(vtype=GRB.BINARY,
                                   name='X_' + str(p) + '_' + str(s) + '_' + str(h) + '_' + str(d) + '_' + str(t))
         for s in S for h in H for d in D for p in Lambda_sd[s, d] for t in Times
         if d <= Theta_p[p] and t + T_ps[p, s] <= nTime  # Before due date, enough time in day to complete surgery
         and A_sd[s, d] >= E_ps[p, s]  # Surgeon has enough available time
         and (t == 0 or t >= minDuration)  # No point in starting in this window
         }

    # Number of unused rooms in hospital h, day d, period t-1
    Y = {(h, d, t): m.addVar(vtype=GRB.INTEGER, lb=0.0) for h in H for d in D for t in Times}

    # 1 if surgeon s is in hospital h on day d
    Z = {(s, h, d): m.addVar(vtype=GRB.BINARY) for s in S for h in H for d in D}

    # 1 if room r in hospital h on day d is used, 0 otherwise
    Ron = {(h, d, r): m.addVar(vtype=GRB.BINARY) for h in H for d in D for r in R_h[h]}

    # 1 if room r in hospital-day (h, d) finishes at time t,
    # t is between end of ordinary hours and end of overtime hours
    Roff = {(h, d, t, r): m.addVar(vtype=GRB.BINARY) for h in H for d in D for t in Times[DaySlots - 1:] for r in
            R_h[h]}

    # Speed-up for conservation of flow constraint summations
    # List of operations that start/finish at time t in (h, d)
    OpStart = {(h, d, t): [] for h in H for d in D for t in Times}
    OpFinish = {(h, d, t): [] for h in H for d in D for t in Times}

    for p, s, h, d, t in X:
        OpStart[h, d, t].append(X[p, s, h, d, t])
        OpFinish[h, d, t + T_ps[p, s]].append(X[p, s, h, d, t])

    del p, s, h, d, t

    # Speed-up for objective function - optional patient summation, patients operated on if due date is after horizon
    OptionalPatients = {(p, s, h, d, t): X[p, s, h, d, t] for p, s, h, d, t in X if Theta_p[p] > max(D)}
    PSHD_List = collections.defaultdict(list)
    for k in X:
        PSHD_List[k[:-1]].append(X[k])

    # --------------------------- OBJECTIVE ----------------------------------------
    m.setObjective(
        quicksum(U_p[p] * X[p, s, h, d, t] for p, s, h, d, t in X)  # Value of scheduling operations
        - quicksum(
            quicksum(
                Ron[h, d, r] * K_hdr[h, d, r]  # Cost of using rooms
                + C_hdr[h, d, r] / (60 / t_period) * quicksum(
                    Roff[h, d, t, r] * (t - DaySlots) for t in Times[DaySlots:])  # Cost of overtime
                for r in R_h[h])
            + quicksum(L_shd[s, h, d] * Z[s, h, d] for s in S)  # Surgeon costs
            for h in H for d in D),
        GRB.MAXIMIZE)

    # --------------------------- CONSTRAINTS ----------------------------------------

    ConservationOfORFlow = {(h, d, t):
        m.addConstr(
            Y[h, d, t] == Y[h, d, t - 1] -  # No of rooms unused = no. in prev time period
            quicksum(x for x in OpStart[h, d, t - 1]) +  # Minus operations starting in t-1
            quicksum(x for x in OpFinish[h, d, t - 1]) -  # Plus operations finishing in t-1
            # Minus rooms finishing for day in this time period
            (quicksum(Roff[h, d, t - 1, r] for r in R_h[h]) if t >= DaySlots + 1 else 0))
        for h in H for d in D for t in Times[1:]}

    # These two constraints enforce mandatory surgeries being performed,
    # and optional surgeries being performed at most once. This also ensures patients are in at most one hospital.
    Mandatory = {p: m.addConstr(quicksum(X[k] for k in X if k[0] == p) == 1)
                 for p in P if Theta_p[p] <= max(D)}

    # Prioritise cheaper rooms
    OrderHDR = {(h, d, r, rd):
                    m.addConstr(Ron[h, d, r] <= Ron[h, d, rd])
                for h in H for d in D for r in R_h[h] for rd in R_h[h]
                if rd != r and (K_hdr[h, d, rd] < K_hdr[h, d, r] or
                                (K_hdr[h, d, rd] == K_hdr[h, d, r] and rd < r))
                }

    # Set up branching variables
    # Optional OP on or not
    W = {p: m.addVar(vtype=GRB.BINARY) for p in P if Theta_p[p] > max(D)}
    Optional = {p: m.addConstr(quicksum(X[k] for k in X if k[0] == p) == W[p]) for p in P if Theta_p[p] > max(D)}

    # Surgeon used on a day or not
    WS = {(s, d): m.addVar(vtype=GRB.BINARY) for s in S for d in D}

    # Dictionary of variables for patient on hospital-day
    PHD_List = collections.defaultdict(list)
    for k in X:
        PHD_List[k[0], k[2], k[3]].append(X[k])
    # Patient on a hospital day
    WP = {k: m.addVar(vtype=GRB.BINARY) for k in PHD_List}
    PHDCon = {k: m.addConstr(WP[k] == quicksum(PHD_List[k])) for k in PHD_List}

    # Number of rooms for each hospital day
    WC = {(h, d): m.addVar(vtype=GRB.INTEGER) for h in H for d in D}
    SetWC = {(h, d): m.addConstr(WC[h, d] == quicksum(Ron[h, d, r] for r in R_h[h]))
             for h in H for d in D}
    for k in W:
        W[k].BranchPriority = 100
    for k in WS:
        WS[k].BranchPriority = 30
    for k in WP:
        WP[k].BranchPriority = 20
    for k in WC:
        WC[k].BranchPriority = 15
    # # for k in Z:
    #     Z[k].BranchPriority = 40
    # for k in Ron:
    #     Ron[k].BranchPriority = 20

    # Stop Gurobi from presolving out these constraints
    m.setParam('GURO_PAR_MINBPFORBID', 1)

    # Surgeons cannot perform more surgeries than they have time
    SurgeonLimit = {(s, h, d):
        m.addConstr(
            quicksum(E_ps[k[0], s] * X[k] for k in X if k[1:-1] == (s, h, d)) <= A_sd[s, d] * Z[s, h, d])
        for s in S for h in H for d in D if A_sd[s, d] > 0}

    # Dictionaries of surgeons, days, and start times, and dict of surgeon, day, end time
    SurgeonStart = collections.defaultdict(list)
    SurgeonEnd = collections.defaultdict(list)
    for p, s, h, d, t in X:
        SurgeonStart[s, d, t].append(X[p, s, h, d, t])
        SurgeonEnd[s, d, t + E_ps[p, s]].append(X[p, s, h, d, t])

    # Availability of surgeon s on day d at time t
    SW = {(s, d, t): m.addVar() for s, d, t in SurgeonStart if t > 0}
    # Make sure surgeon is only doing one surgery at a time
    CountSurgeon = {(s, d, t):
                        m.addConstr(SW[s, d, t] == (SW[s, d, t - 1] if t > minDuration else
                                                    quicksum(Z[s, h, d] for h in H) - quicksum(SurgeonStart[s, d, 0]))
                                    + quicksum(SurgeonEnd[s, d, t]) - quicksum(SurgeonStart[s, d, t]))
                    for s, d, t in SW}

    # Surgeons can only be in one hospital per day, and only if they are active that day
    SurgeonInOneHospital = {(s, d): m.addConstr(quicksum(Z[s, h, d] for h in H) == WS[s, d]) for s in S for d in D}

    # This enforces above constraints
    # Surgeon can't operate in hospital h unless they are at that hospital
    LinkXToZ = {k: m.addConstr(quicksum(PSHD_List[k]) <= Z[k[1:]]) for k in PSHD_List}

    # Initial number of unused ORs on (h, d) = number of rooms turned on
    InitialY = {(h, d): m.addConstr(Y[h, d, 1] == quicksum(Ron[h, d, r] for r in R_h[h])) for h in H for d in D}

    # # If a room is turned on, it must be turned off
    ConservationOfOnOff = {(h, d, r): m.addConstr(
        Ron[h, d, r] == quicksum(Roff[h, d, t, r] for t in Times[DaySlots - 1:])
    ) for h in H for d in D for r in R_h[h]}

    m.optimize()

    # --------------------------- RESULTS ----------------------------------------
    with open('V' + version + '_' + str(t_period) + 'min_' + folder_name + '.txt', 'a') as f:

        for d in D:
            string = '\nDay ' + str(d) + ':\n'
            f.write(string)
            for h in H:
                # If any operations on (h, d), check if Ron variables work
                if any([X[k].x for k in X if k[2:-1] == (h, d)]):
                    # print hospital & rooms available
                    # print('Hospital ', h, ', rooms: ', [round(abs(Y[h, d, t].x), 2) for t in Times])
                    string = 'Hospital ' + str(h) + '\n'
                    f.write(string)
                    for r in R_h[h]:
                        string = 'room ' + str(r) + ' in use: ' + str(round(abs(Ron[h, d, r].x), 2)) + '. Time off: ' + \
                                 str([t for t in Times[DaySlots - 1:] if Roff[h, d, t, r].x > 0.9]) + '\n'
                        f.write(string)
                    # Explain what operation info is given
                    f.write('(p, s, [t_start, t_end])\n')
                # Print operations scheduled on (h, d) as (patient, surgeon, start time slot, end time slot)
                for t in Times:
                    for p in P:
                        for s in S:
                            if (p, s, h, d, t) in X and X[p, s, h, d, t].x > 0.9:
                                string = str((p, s, (t, t + T_ps[p, s]))) + '\n'
                                f.write(string)

    # --------------------------- ERROR CHECKING ----------------------------------------

    # Mandatory patients who are not scheduled
    missed = [p for p in P if not (any([X[k].x for k in X])) and Theta_p[p] <= max(D)]
    if len(missed) > 1:
        print('ERROR - Missed mandatory patients = ', missed)

    # Surgeons with no free time and do operations
    nosurgeontime = [(s, d) for s in S for d in D if (any([X[k].x for k in X
                                                           if k[1] == s and k[-2] == d])) and A_sd[s, d] == 0]

    if len(nosurgeontime) > 0.1:
        print('ERROR - Surgeons working with no time = ', nosurgeontime)

    # Surgeons operating in multiple hospitals
    SurgeonMultipleHospital = {(s, d): [h for h in H if any([Z[s, h, d].x > 0.1])] for s in S for d in D}

    for k in SurgeonMultipleHospital:
        if len(SurgeonMultipleHospital[k]) > 1.1:
            print('ERROR - Surgeon s working multiple hospitals on day d = ', k)

    # Rooms finished more than once
    RoomMultipleFinish = {(h, d, r): [t for t in Times[DaySlots - 1:] if any([Roff[h, d, t, r].x > 0.1])] for d in D for
                          h
                          in H for r in R_h[h]}

    for h, d, r in RoomMultipleFinish:
        if len(RoomMultipleFinish[h, d, r]) > 1.1:
            print('ERROR - Room ', r, ' finishing multiple times on hospital-day ', (h, d))

    return


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# Optimisation

# Set up command line arguments
# User can enter up to four - no. of patients, trial no., time period length, and absolute path of CORPS-Data
# Number of patients
size = 35
if len(sys.argv) > 1.1:
    size = int(sys.argv[1])
    if size not in [35, 40, 45, 70, 80, 90, 105, 120, 135]:
        raise ValueError('Not a valid patient set size.')

# Trial number
trial = 1
if len(sys.argv) > 2.1:
    trial = int(sys.argv[2])
    if trial not in range(1, 6):
        raise ValueError('Not a valid trial number.')

# Length of time periods
t_step = 10
if len(sys.argv) > 3.1:
    t_step = float(sys.argv[3])
    if t_step <= 0:
        raise ValueError('Not a valid time period length.')

# Absolute path to CORPS-Data folder (e.g. C:\Users\admin\OR')
data_loc = r'C:\Users\admin\Documents\University\OR_Project'
if len(sys.argv) > 4.1:
    data_loc = sys.argv[4]

# Optimise
surgery_optimise(120, trial, t_step, data_loc)