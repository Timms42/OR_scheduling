# Author: Liam Timms
# Change log: Version 4
# Prev. version - allocation_reader_V3.py
# Automatically detect the room for each scheduled operation.

from re import findall

def read_allocation(file_name):
    """
    :param file_name: (str) Absolute path & file name of the output allocation text file.
    :return: (dict) Dict of dicts, and dict.
                    First contains of patient-surgeon-time allocations for each hospital and day.
                    Indexed by dict[hospital][day]
                    Second contains room usages for each hospital and day.
                    Indexed by dict[hospital][day]
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
    hd_schedule = {(h, d): [] for h in range(1, nhosp+1) for d in range(1, ndays+1)}

    # Store h,d room usage in dictionary (Note: h, d start at 1 here)
    hd_rooms = {(h, d): [] for h in range(1, nhosp + 1) for d in range(1, ndays + 1)}

    # For each day in the horizon (labelled 0,..., ndays-1 for ease of indexing)
    for d in range(ndays):
        # Subset of allocation file corresponding to day d
        # If 1st day, slice from start to beginning of day 2
        if d == 0:
            lines_d = lines[:d_ind[d+1]]

        # If last day, slice from end of 2nd last day up to end
        elif d == ndays-1:
            lines_d = lines[d_ind[-1]:]

        # Otherwise, slice from day d up to beginning of day d+1
        else:
            lines_d = lines[d_ind[d]:d_ind[d+1]]

        # Get indices of all 'Hospital h' lines in allocation file
        # Index is 'None' if hospital has no operations on that day
        h_ind = [lines_d.index(hospitals[h]) if hospitals[h] in lines_d else 'None' for h in range(len(hospitals))]

        # Index over each hospital, labelled 0, 1, ..., n for ease of indexing
        for h in range(nhosp):
            # If the hospital has operations on this day, get (h, d) subsection of allocation text
            if h_ind[h] != 'None':
                # If all later hospitals have no operations, take slice of Lines_d from 'Hospital h:' to the end of file
                if all([hospital == 'None' for hospital in h_ind[h+1:]]):
                    lines_hd = lines_d[h_ind[h]:]

                # If the next hospital has no operations, find the next hospital to have operations
                elif h_ind[h+1] == 'None':
                    # Get index of operations
                    next_hosp = [ii for ii in h_ind[h+1:] if ii != 'None'][0]
                    lines_hd = lines_d[h_ind[h]:next_hosp]

                # Otherwise, take slice of lines_d from 'Hospital h' to 'Hospital h+1'
                else:
                    lines_hd = lines_d[h_ind[h]:h_ind[h+1]]

            # If the hospital has no operations on this day, the (h, d) subsection text is empty
            elif h_ind[h] == 'None':
                lines_hd = []

            else:
                raise ValueError('Something is wrong in checking the hospital data for (h,d) ', h, d)   # ERROR CHECKING

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

                hd_rooms[h + 1, d + 1] = rooms_data

                # Remove irrelevant information by slicing lines_hd after '(p, s, [t_start, t_end])'
                lines_hd = lines_hd[alloc_key+1:]

                # Remove any blank lines read in
                if '' in lines_hd:
                    lines_hd.remove('')

                # Split strings into list of [p, s, t_start, t_start+T_ps]
                for o in range(len(lines_hd)):
                    lines_hd[o] = lines_hd[o][1:-1].split(',')

                    for ind in range(len(lines_hd[o])):
                        # Remove surrounding brackets and whitespace
                        lines_hd[o][ind] = float(lines_hd[o][ind].replace(' ', '').replace('(', '').replace(')', ''))

            # Save hd subsection to dictionary
            hd_schedule[h+1, d+1] = lines_hd

    # Make empty dicts to store individual hospital allocations & room usages
    indiv_alloc = {hs: 0 for hs in range(1, nhosp+1)}
    indiv_rooms = {hs: 0 for hs in range(1, nhosp + 1)}

    for hkey in indiv_alloc:
        # Make dictionary of hospital 'hkey' allocations over the horizon
        # Indexed indiv_alloc[hospital][day], returns list of operation allocations
        indiv_alloc[hkey] = {k[1]: hd_schedule[k] for k in hd_schedule.keys() if k[0] == hkey}
        indiv_rooms[hkey] = {k[1]: hd_rooms[k] for k in hd_rooms.keys() if k[0] == hkey}

    return indiv_alloc, indiv_rooms


read_ops, read_rooms = read_allocation(r'C:\Users\admin\Documents\University\2021\COSC3000\surgery_allocation.txt')

ops = read_ops[2][3]
rooms = read_rooms[2][3]

ntime = 90
# List of 1 if room r is available in time period t, 0 otherwise
available = {r[0]: [1]*ntime for r in rooms}


# TAKE 2 - generalising room list
for ii in range(len(ops)):
    # Start and end time periods of operation ops[ii]
    tstart = int(ops[ii][2])
    tend = int(ops[ii][3])

    for r in rooms:
        # Get the room number
        rnum = r[0]
        # If room r is available for the whole operation in AND operation hasn't been assigned a room
        # (i.e. it only has 4 elements), append r to the end of the operation info.
        # Note: end of op ii and start of op jj can be the same, since we assume end is at the start of the time period.
        #       So, we don't need to check if available at tstart
        if all(available[rnum][tstart+1:tend+1]) and len(ops[ii]) < 4.1:
            ops[ii].append(rnum)
            for jj in range(tstart, tend+1):
                available[rnum][jj] = 0

print(ops)