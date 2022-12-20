import pylab as pyl
import numpy as np
import matplotlib as mp
import matplotlib.cm as cm
import matplotlib.colors as clr
import pandas as pd
from re import findall
from plotly.express import timeline
from collections import defaultdict
from PIL import ImageColor

# Change log: Version 15
# Modified: 6 Jan 2022
# Tidying up figures
# Make row labels neater

# Reference: Colors from www.colorbrewer2.org by Cynthia A. Brewer, Geography, Pennsylvania State University.

# ----------------------------------------------------------
# Timms and Forbes (2021) UQ operations research project
# Import data
# Data from Roshaenaei et al. (2018)


def read_file(absolute, folder, name):
    # abs (string): absolute location of folder
    # folder (string): name of folder with data set
    # name (string): file name of data set without extension , e.g. 'A_sd'

    # Create list of variables in the file
    index = []
    for v in name.partition('_')[2]:
        index.append(v)

    # Read in csv file as dataframe
    df = pd.read_csv(absolute + "\\" + folder + "\\" + name + '.csv', index_col=index)

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


def read_allocation(file_name, num_periods, room_costs):
    """
    :param file_name: (str) Absolute path & file name of the output allocation text file.
    :param num_periods: (int) Number of time periods
    :return:    (dict) Dict of dicts. Contains of patient-surgeon-time allocations for each hospital and day.
                    Indexed by dict[hospital][day]
                (dict) Room usages for each hospital and day.
                    Indexed by dict[hospital][day]
                (int) Number of days in file_name
    """
    with open(file_name, 'r') as f:
        # Get list of lines in text file
        lines = f.readlines()
        # Strip whitespace and formatting characters
        lines = [l.strip() for l in lines]

    # Find all 'Day d:' strings in text file, and isolate the day numbers
    days = [s for s in lines if 'Day' in s]

    # Set number of days
    ndays = len(days)

    # Indices of all 'Day d:' lines in file
    d_ind = [lines.index(d) for d in days]

    # Find all 'Hospital h:' strings in text file, and isolate the hospital numbers.
    hospitals = list(set([s for s in lines if 'Hospital' in s]))
    # Sort in ascending order
    hospitals.sort()
    # Set number of days
    nhosp = len(hospitals)

    # Store h,d subsections in dictionary (Note: h, d start at 1 here)
    hd_schedule = {(h, d): [] for h in range(1, nhosp + 1) for d in range(1, ndays + 1)}

    # Store h,d room usage in dictionary (Note: h, d start at 1 here)
    hd_rooms = {(h, d): [] for h in range(1, nhosp + 1) for d in range(1, ndays + 1)}

    # For each day in the horizon (labelled 0,..., ndays-1 for ease of indexing)
    for d in range(ndays):
        # Subset of allocation file corresponding to day d
        # If 1st day, slice from start to beginning of day 2
        if d == 0:
            lines_d = lines[:d_ind[d + 1]]

        # If last day, slice from end of 2nd last day up to end
        elif d == ndays - 1:
            lines_d = lines[d_ind[-1]:]

        # Otherwise, slice from day d up to beginning of day d+1
        else:
            lines_d = lines[d_ind[d]:d_ind[d + 1]]

        # Get indices of all 'Hospital h' lines in allocation file
        # Index is 'None' if hospital has no operations on that day
        h_ind = [lines_d.index(hospitals[h]) if hospitals[h] in lines_d else 'None' for h in range(len(hospitals))]

        # Index over each hospital, labelled 0, 1, ..., n for ease of indexing
        for h in range(nhosp):
            # If the hospital has operations on this day, get (h, d) subsection of allocation text
            if h_ind[h] != 'None':
                # If all later hospitals have no operations, take slice of Lines_d from 'Hospital h:' to the end of file
                if all([hospital == 'None' for hospital in h_ind[h + 1:]]):
                    lines_hd = lines_d[h_ind[h]:]

                # If the next hospital has no operations, find the next hospital to have operations
                elif h_ind[h + 1] == 'None':
                    # Get index of operations
                    next_hosp = [ii for ii in h_ind[h + 1:] if ii != 'None'][0]
                    lines_hd = lines_d[h_ind[h]:next_hosp]

                # Otherwise, take slice of lines_d from 'Hospital h' to 'Hospital h+1'
                else:
                    lines_hd = lines_d[h_ind[h]:h_ind[h + 1]]

            # If the hospital has no operations on this day, the (h, d) subsection text is empty
            elif h_ind[h] == 'None':
                lines_hd = []

            else:
                raise ValueError('Something is wrong in checking the hospital data for (h,d) ', h, d)  # ERROR CHECKING

            # If the hd subsection is not an empty list, extract patient-surgeon allocations and times
            if lines_hd:
                # Find the location of the allocation key in the hd subsection
                alloc_key = lines_hd.index('(p, s, [t_start, t_end])')

                # Save room usage data, exclude first line ('Hospital h')
                rooms_data = lines_hd[1:alloc_key]

                # Extract room number, active status, and time off
                for r in range(len(rooms_data)):
                    rooms_data[r] = [findall(r'\d\.\d+|\d+', s) for s in rooms_data[r].split()]

                    rooms_data[r] = [float(num[0]) for num in rooms_data[r] if num]

                # Filter out all inactive rooms, i.e. rooms_data[r][1] = 0
                rooms_data = [item for item in rooms_data if item[1] > 0.1]

                # Sort list of active rooms by ascending cost for that hospital-day (K_hdr indexes from 1)
                rooms_data = sorted(rooms_data, key=lambda x: room_costs[h + 1, d + 1, x[0]])

                hd_rooms[h + 1, d + 1] = rooms_data

                # Remove irrelevant information by slicing lines_hd after '(p, s, [t_start, t_end])'
                lines_hd = lines_hd[alloc_key + 1:]

                # Remove any blank lines read in
                if '' in lines_hd:
                    lines_hd.remove('')

                # Split strings into list of [p, s, t_start, t_start+T_ps]
                for o in range(len(lines_hd)):
                    lines_hd[o] = lines_hd[o][1:-1].split(',')

                    for ind in range(len(lines_hd[o])):
                        # Remove surrounding brackets and whitespace
                        lines_hd[o][ind] = float(lines_hd[o][ind].replace(' ', '').replace('(', '').replace(')', ''))

                # List of 1 if room r is available in time period t, 0 otherwise
                available = {r[0]: [1] * num_periods for r in rooms_data}

                for ii in range(len(lines_hd)):
                    # Start and end time periods of operation ops[ii]
                    tstart = int(lines_hd[ii][2])
                    tend = int(lines_hd[ii][3])

                    for r in rooms_data:
                        # Get the room number
                        rnum = int(r[0])

                        # If room r is available for the whole operation in AND operation hasn't been assigned a room
                        # (i.e. it only has 4 elements), append r to the end of the operation info.
                        # Note: end of op ii and start of op jj can be the same, since we assume end
                        # is at the start of the time period. So, we don't need to check if available at tstart
                        if all(available[rnum][tstart + 1:tend + 1]) and len(lines_hd[ii]) < 4.1:
                            lines_hd[ii].append(rnum)
                            for jj in range(tstart, tend):
                                available[rnum][jj] = 0

            # Save hd subsection to dictionary
            hd_schedule[h + 1, d + 1] = lines_hd

    # Make empty dicts to store individual hospital allocations & room usages
    indiv_alloc = {hs: 0 for hs in range(1, nhosp + 1)}
    indiv_rooms = {hs: 0 for hs in range(1, nhosp + 1)}

    for hkey in indiv_alloc:
        # Make dictionary of hospital 'hkey' allocations over the horizon
        # Indexed indiv_alloc[hospital][day], returns list of operation allocations
        indiv_alloc[hkey] = {k[1]: hd_schedule[k] for k in hd_schedule.keys() if k[0] == hkey}
        indiv_rooms[hkey] = {k[1]: hd_rooms[k] for k in hd_rooms.keys() if k[0] == hkey}

    return indiv_alloc, indiv_rooms, ndays


def ConvertColourInput(cList):
    """
    :param cList:  EITHER (str) name of colormap from matplotlib cm
                    OR (list) list of RGB triplets or HEX codes
    :return:    (list) list of colours as RGBA tuples
    """
    # If the 'colour' input is a colormap name or single HEX code
    if isinstance(cList, str):
        # Check if single HEX code, i.e. begins with #
        if cList[0] == '#':
            return cList
        else:
            # Unpack colormap and take 3 colours to assign to each different room
            colourmap = cm.get_cmap(cList, 100)
            colourlist = [colourmap(x) for x in [30, 40, 50, 60, 65, 75, 80, 85, 90, 95]]

    # Otherwise, if the 'colour' input is a list
    elif isinstance(cList, list):
        # If the elements are RGB triplets
        if all([isinstance(c, float) for c in cList]):
            # Check if RGB triplets were given as whole numbers or between [0, 1]
            if not all([-0.1 < num < 1.1 for RGB in cList for num in RGB]):
                # Convert to RGB values to between 0 and 1
                colourlist = [[num / 255 for num in RGB] for RGB in cList]
            else:
                colourlist = cList
        # If the elements are HEX codes, convert to RGBA
        elif all([isinstance(c, str) for c in cList]):
            # colourlist = [clr.to_rgba(c, 1) for c in cList]
            colourlist = cList
        else:
            raise ValueError('Input "colour" is not a colormap name nor a list of RGB/HEX triplets.')

    else:
        raise ValueError('Input "colour" is not a colormap name nor a list of RGB/HEX triplets.')

    return colourlist


def CreateGantt(alloc, cList, transparency, hospitals, d, d_names, duration_data, prep_data, clean_data,
                t_step, numT, rList, sList, numP):
    """
    :param alloc: (dict) operation allocations for all hospitals and days
    :param cList: EITHER (str) name of colormap from matplotlib cm
                    OR (list) list of RGB triplets
    :param transparency: (float) transparency value between 0 and 1 for colour of prep/clean time bars
    :param hospitals: (list) list of hospital to plot
    :param d: (int) day number
    :param d_names: (list) list of the names of days (as strings)
    :param duration_data: (dict) dictionary of surgery durations, indexed by [patient, surgeon]
    :param prep_data: (dict) dictionary of surgery prep durations, indexed by [patient]
    :param clean_data: (dict) dictionary of post-op cleaning durations, indexed by [patient]
    :param cList: (list)
    :param t_step: (float) length of time periods
    :param numT: (int) total number of time periods
    :param rList: (list) list of the rooms in the hospital
    :param sList: (list) list of the surgeons in the hospital
    :param numP: (int) number of patients to be scheduled

    :return: figure object to plot
    """

    # Calculate start and end of the day as a timestamp
    morning = 8  # 8am start of the day
    sod = pd.to_datetime(morning, unit='h') # Start of day time
    morning_offset = pd.Timedelta(morning, unit='h')  # Add this to all times - times are relative to start of day
    eod = pd.to_datetime(numT * t_step, unit='m') + morning_offset  # End of day time

    # Convert RGB tuples from in [0, 1] to in [0, 255], keep transparency between [0, 1]
    # cList = [tuple([c * 255 for c in rgb_tup[:-1]] + [rgb_tup[-1]]) for rgb_tup in cList]

    # Create dictionary linking surgeons to colours
    colourDict = dict(zip(sList, cList))

    # Initialise list for arc dictionaries
    timeline_list = []

    # Intialise list for opacities
    opacity_list = []

    # Get the (h, d) subset
    for h in hospitals:
        ops = alloc[h][d]

        for arc in ops:
            # Check if list is not empty
            if arc:
                # Take p, s, t_start from the read allocation data
                p, s, t_start, r = int(arc[0]), int(arc[1]), arc[2], int(arc[4])
                # Convert relative (to the start of day) times to absolute start/end times
                # Start times of preparation, operation, cleaning, and end time
                prep_start = t_start * t_step

                op_start = prep_start + prep_data[p]

                clean_start = op_start + duration_data[p, s]

                t_end = clean_start + clean_data[p]

                prep_start = pd.to_datetime(prep_start, unit='m') + morning_offset
                op_start = pd.to_datetime(op_start, unit='m') + morning_offset
                clean_start = pd.to_datetime(clean_start, unit='m') + morning_offset
                t_end = pd.to_datetime(t_end, unit='m') + morning_offset

                # Create room, surgeon, and hospital labels for this arc
                room = 'Room ' + str(r)
                surg = 'Surgeon ' + str(s)
                hosp = str(h)

                # prep_colour = tuple([num for num in cList[r-1][:-1]] + [transparency])
                # Convert RGBA tuples to HEX codes for timeline
                # prep_colour = clr.to_hex(prep_colour, True)
                # op_colour = clr.to_hex(cList[r-1], True)

                # op_colour = 'rgba' + str(cList[r-1])
                # prep_colour = 'rgba' + str(prep_colour)

                # Prep goes from start to start of operation. Op goes from start of op to start of clean.
                # Clean goes from end of op to end of job.
                prep_dict = defaultdict(patient=p, Surgeon=surg, start=prep_start, finish=op_start, Rooms=room,
                                        Hospital=hosp, colour=cList[r-1], delta=prep_data[p],
                                        category='Preparation')
                mainop_dict = defaultdict(patient=p, Surgeon=surg, start=op_start, finish=clean_start, Rooms=room,
                                          Hospital=hosp, colour=cList[r-1], delta=duration_data[p, s],
                                          category='Operation')
                clean_dict = defaultdict(patient=p, Surgeon=surg, start=clean_start, finish=t_end, Rooms=room,
                                         Hospital=hosp, colour=cList[r-1], delta=clean_data[p],
                                         category='Cleaning')

                # Update list of operations in timeline, add corresponding transparencies
                for job_dict in [prep_dict, mainop_dict, clean_dict]:
                    timeline_list.append(job_dict)
                    opacity_list += [transparency, 1, transparency]

    gantt_df = pd.DataFrame(timeline_list)
    # Create multiline title
    the_title = 'CORPS Schedule: ' + d_names[d] + '<br> Total patients:' + str(numP) +\
                '. Total surgeons: ' + str(len(sList))
    # Plot the Gantt chart
    # facet_row creates subplots for each hospital
    # category_orders forces the y_axis tp be sorted according to the given list
    f = timeline(gantt_df, x_start='start', x_end='finish', y='Rooms', range_x=[sod, eod], color='Surgeon',
                 facet_row='Hospital',
                 category_orders={'Rooms': rList, 'Surgeon': sList}, title=the_title,
                 color_discrete_map=colourDict, opacity=opacity_list)

    # Change box border colours to black for contrast
    f.update_traces(marker=dict(line=dict(color='black')))

    return f, gantt_df


def AllocationToCharts(filename, dList, dNames, dur_data, prep_data, clean_data, cost_data, t_step,
                       numT, numP, rList, sList, cList, transparency, backCol, save):
    """
    :param filename: (str) absolute path of allocation data file
    :param dList: list of the days to graph
    :param dNames: dictionary of (day number: day names) pairs
    :param dur_data: (dict) dictionary of surgery durations, indexed by [patient, surgeon]
    :param prep_data: (dict) dictionary of surgery prep durations, indexed by [patient]
    :param clean_data: (dict) dictionary of post-op cleaning durations, indexed by [patient]
    :param cost_data: (dict) dictionary of room cost/hour, indexed by [hospital, day, room]
    :param t_step: (float) length of time periods
    :param numT: (int) total number of time periods
    :param numP: (int) number of patients
    :param rList: rList: (list) list of the rooms in the hospital
    :param sList: (list) list of the surgeons in the hospital
    :param cList: : EITHER (str) name of colormap from matplotlib cm
                    OR (list) list of RGB triplets / HEX codes. List must have 10 colours.
    :param transparency: (float) transparency value between 0 and 1 for colour of prep/clean time bars
    :param backCol: EITHER (str) name of colour from matplotlib cm
                    OR (list) RGB triplet
                    OR (str) HEX code
    :param save: (int) 1 to save chart as png, 0 to not save
    :return: None, create Gantt charts for all days in dList based on alloc
    """

    # Read in the data
    allocate_ops, allocate_rooms, nDays = \
        read_allocation(filename, numT, cost_data)

    # Convert colour list and background colour to appropriate colormap or RGB/Hex list
    cList = ConvertColourInput(cList)
    backCol = ConvertColourInput(backCol)

    for input_day in dList:
        # List of hospitals with operations on that hospital-day
        active_hospitals = [h for h in [1, 2, 3] if allocate_ops[h][input_day]]

        # -----------------------------PLOTTING THE GANTT CHART--------------------------------------
        # Create and show the Gantt chart for day input_day
        fig, the_df = CreateGantt(allocate_ops, cList, transparency, active_hospitals, input_day, dNames, dur_data,
                                  prep_data, clean_data, t_step, numT, rList, sList, numP)

        # Change background colour of the graph
        fig.update_layout(plot_bgcolor=backCol, xaxis_title='Time')

        # Set xticks at every hour
        fig.update_xaxes(dtick=60 * 60 * 1000,
                         tickformat='%H:%M')

        savename = 'gantt_chart_size' + str(numP) + '_day' + str(input_day)

        if save == 'png':
            # .png friendly font size
            fig.update_layout(font=dict(size=25))
            # Save graph as .png
            savename += '.png'
            fig.write_image(savename, width=1980, height=1080)

        elif save == 'html':
            # Webpage friendly font size
            fig.update_layout(font=dict(size=16))
            savename += '.html'
            fig.write_html(savename, full_html=False, include_plotlyjs='cdn')

        if save == 'pdf':
            # .pdf friendly font size
            fig.update_layout(font=dict(size=25))
            # Save graph as .pdf
            savename += '.pdf'
            fig.write_image(savename, width=1980, height=1080)

        else:
            # Do not save the figure
            fig.update_layout(font=dict(size=16))

        fig.show()

    return the_df


# -----------------------------IMPORT THE DATA----------------------------------
# Name of subfolder with data set
patient_no = 120
trial_no = 1

file_loc = r'C:\Users\admin\Documents\University\OR_Project\CORPS-Data\Alpha'
folder_name = 'size_' + str(patient_no) + '_trial_' + str(trial_no)

# Read in data files
A_sd = read_file(file_loc, folder_name, 'A_sd')
B_hdr = read_file(file_loc, folder_name, 'B_hdr')
C_hdr = read_file(file_loc, folder_name, 'C_hdr')  # Note: C_hdr = 3/16 * K_hdr
Delta_sd = read_file(file_loc, folder_name, 'Delta_sd')
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
# Round times to nearest 10 minutes
t_period = 10
# Number of time periods
nTime = 60
Tmins_ps = T_ps.copy()  # Preserve a copy of T_ps in minutes

for dur in E_ps:
    E_ps[dur] = dur2period(E_ps[dur], t_period) * 10  # Actual surgery duration
    T_ps[dur] = dur2period(T_ps[dur], t_period) * 10  # Total sugery duration = F_p + E_ps + G_p

# Round prep and cleaning time F_p, G_p
for p_index in F_p:
    F_p[p_index] = dur2period(F_p[p_index], t_period) * 10  # Prep time
    G_p[p_index] = dur2period(G_p[p_index], t_period) * 10  # Clean time

# -----------------------------THE PROGRAM--------------------------------------
# Choose colormap and list of colours
colormap = 'Blues'
colours = ['#feebe2', '#fbb4b9', '#f768a1', '#c51b8a', '#7a0177', '#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026']

# Create list of days, rooms, and surgeons (assuming rooms and surgeons are homogenous across hospitals and days
days = [d for d in range(1, max([i[1] for i in A_sd.keys()]) + 1)]

# TESTING - use 1st day
days = [1]

day_names = {1: 'Monday',
             2: 'Tuesday',
             3: 'Wednesday',
             4: 'Thursday',
             5: 'Friday',
             6: 'Saturday',
             7: 'Sunday'}

rooms = ['Room ' + str(r) for r in range(1, max([i[2] for i in K_hdr.keys()]) + 1)]

surgeons = ['Surgeon ' + str(s) for s in range(1, max([i[0] for i in A_sd.keys()]) + 1)]

datafile = r'C:\Users\admin\Documents\University\2021\COSC3000\Visualisation_project\surgery_allocation_' \
           + folder_name + '.txt'

# Transparency for prep/clean time bars
prep = 1/2


# Parameter 'save' must be string, either 'png', 'pdf', 'html'
the_dataframe = AllocationToCharts(datafile, days, day_names, E_ps, F_p, G_p, K_hdr, t_period, nTime, patient_no,
                                   rooms, surgeons, colours, prep, '#677E8C', save='pdf')
