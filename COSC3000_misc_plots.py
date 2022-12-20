import pylab
import numpy as np
import cmasher as cmr
import matplotlib.pyplot as plt
from matplotlib import cm
from pandas import read_csv, DataFrame
from collections import defaultdict
from plotly.express import scatter, histogram
import plotly.figure_factory as ff


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


# Name of subfolder with data set
file_loc = r'C:\Users\admin\Documents\University\OR_Project\CORPS-Data\Alpha'
folder_name = 'size_35_trial_4'

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
nTime = 90
Tmins_ps = T_ps.copy()  # Preserve a copy of T_ps in minutes

for dur in E_ps:
    E_ps[dur] = dur2period(E_ps[dur], t_period)
    T_ps[dur] = dur2period(T_ps[dur], t_period)

# Round prep time F_p
for p_index in F_p:
    F_p[p_index] = dur2period(F_p[p_index], t_period)

# --------------------------THE PROGRAM---------------------------------
# Generate relevant sets
days = [d for d in range(1, max([i[1] for i in A_sd.keys()]) + 1)]

hospitals = [h for h in range(1, max([i[0] for i in K_hdr.keys()]) + 1)]

rooms = [r for r in range(1, max([i[2] for i in K_hdr.keys()]) + 1)]

surgeons = [s for s in range(1, max([i[0] for i in A_sd.keys()]) + 1)]

patients = [p for p in range(1, max([i[0] for i in E_ps.keys()]) + 1)]

# Colour list
colours = ['#feebe2', '#fbb4b9', '#f768a1', '#c51b8a', '#7a0177', '#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026']
colourmap = cmr.ember

# Create list of surgeon-colour pairs
colourDict = dict(zip(surgeons, colours))

# Data set examined
dataset = Tmins_ps
dataname = 'T_ps'

# Convert Tmins_ps to a pandas dataframe, for the scatter plot
scatter_list = []
Tsurgeons = {s: [] for s in surgeons}
for k in dataset:
    duration_dict = defaultdict(Patient=k[0], Surgeon='Surgeon ' + str(k[1]), Duration=dataset[k])
    scatter_list.append(duration_dict)

    Tsurgeons[k[1]].append(defaultdict(patient=k[0], Duration=dataset[k]))

# Convert dictionary into pandas dataframe for scatter/hist graphing
Tps_df = DataFrame(scatter_list)

# Calculate relevant statistics
data_mean = np.mean([dataset[k] for k in dataset])

# Create np.array() for plotting K_hdr room costs
cost_data = {d:
    np.array(
        [
            [K_hdr[h, d, r] for r in rooms]
            for h in hospitals]
    ) for d in days}

# ---------------CREATE PLOTS-------------------------
# ---------------FIGURE 1 - SCATTERPLOT of surgery durations, coloured by surgeons-------------------------
# # Sort the dataframe by longest surgery duration
scatter_df = Tps_df.copy()
scatter_df.sort_values('Duration', ascending=False, inplace=True)
# Change patient names to strings to allow categorical axes and sorting of axis
scatter_df['Patient'] = scatter_df['Patient'].astype(str)
# Create list of patients, sorted by longest surgery duration
pList = []
for p in scatter_df['Patient']:
    if p not in pList:
        pList.append(p)
fig1 = scatter(scatter_df, x='Patient', y='Duration', color='Surgeon', color_discrete_map=colourDict,
               title='Scatter plot of surgery durations',
               category_orders={
                   'Patient': pList,
                   'Surgeon': ['Surgeon ' + str(s) for s in surgeons]},
               labels={
                   'Duration': 'Duration (mins)',
                   'Patient': 'Patient no.'
               })
# Increase font size
fig1.update_layout(font=dict(size=20))
# Customise markers
fig1.update_traces(marker=dict(size=12,
                               line=dict(width=2,
                                         color='black'))
                   )
fig1.show()
#
# # Save graph as .png
# savename1 = 'Figures/' + dataname + '_scatter' + '.html'
# # fig1.write_image(savename1, width=1980, height=1080)
# fig1.write_html(savename1, full_html=False, include_plotlyjs='cdn')

# ---------------FIGURE 2 - HISTOGRAM of surgery durations-------------------------
# fig2 = histogram(Tps_df, x='Duration', title='Histogram of surgery durations',
#                  color_discrete_sequence=['#998ec3'], marginal='box')
# fig2.update_layout(bargap=0.1, plot_bgcolor='#f7f7f7', font=dict(size=20))
# fig2.add_vline(x=data_mean, line_width=3, line_dash="dash", line_color='#f1a340', annotation_text='Mean duration',
#                annotation_position="top right")
#
# fig2.show()
#
# savename2 = 'Figures/' + dataname + '_histogram' + '.png'
# fig2.write_image(savename2, width=1980, height=1080)
# # fig2.write_html(savename2, full_html=False, include_plotlyjs='cdn')

# ---------------FIGURE 3 - Histogram coloured by patient or surgeon-------------------------
# fig3 = histogram(Tps_df, x='Duration', title='Histogram of surgery durations, coloured by patient',
#                  color='Patient')
# fig3.update_layout(bargap=0.1, plot_bgcolor='#f7f7f7')
# fig3.add_vline(x=data_mean, line_width=3, line_dash="dash", line_color='black', annotation_text='Mean duration',
#                annotation_position="top right")
#
# fig3.show()
# savename3 = 'Figures/' + dataname + '_histogram_bad_p.png'
# fig3.write_image(savename3, width=1980, height=1080)

# ---------------FIGURE 4 - HEATMAP of room costs-------------------------
# colour_scheme = cmr.amber
# # Get min and max costs
# minK = min(K_hdr.values())
# maxK = max(K_hdr.values())
# # Round min & max costs to nearest multipl of 1000
# minK = np.floor(minK / 1000) * 1000
# maxK = np.ceil(maxK / 1000) * 1000
#
#
# fig4, ax4 = plt.subplots(2, int(np.ceil(len(days) / 2)))
# im0 = ax4[0, 0].imshow(cost_data[1], cmap=colour_scheme, vmin=minK, vmax=maxK)
# im1 = ax4[0, 1].imshow(cost_data[2], cmap=colour_scheme, vmin=minK, vmax=maxK)
# im2 = ax4[1, 0].imshow(cost_data[3], cmap=colour_scheme, vmin=minK, vmax=maxK)
# im3 = ax4[1, 1].imshow(cost_data[4], cmap=colour_scheme, vmin=minK, vmax=maxK)
# im4 = ax4[1, 2].imshow(cost_data[5], cmap=colour_scheme, vmin=minK, vmax=maxK)
#
# the_colourbar = fig4.colorbar(im1, ax=ax4.ravel().tolist())
# # Add ticks for colourbar
# numTicks = 4    # Number of bar ticks
# distTick = (maxK - minK)/numTicks   # Distance between bar ticks
# # Create list of colourbar tick places, including min and max
# tick = [ii * distTick + minK for ii in range(numTicks+1)]
# the_colourbar.set_ticks(tick)
# the_colourbar.set_ticklabels(['$' + str(int(t)) for t in tick])
#
# # Flip the y axes so hospitals increase upwards
# for axis_row in ax4:
#     for axis in axis_row:
#         axis.invert_yaxis()
#
# # Remove 6th subplot
# fig4.delaxes(ax4[0, 2])
#
# for h_set in ax4:
#     for axis in h_set:
#         # Set x ticks and labels - rooms
#         axis.set_xticks(np.arange(len(rooms)))
#         axis.set_xticklabels(rooms)
#
#         # Set y ticks, no labels - hospitals
#         axis.set_yticks(np.arange(len(hospitals)))
#         axis.set_yticklabels([])
#
# # Set the y axis labels for 1st column
# ax4[0, 0].set_yticklabels(hospitals)
# ax4[1, 0].set_yticklabels(hospitals)
# ax4[0, 0].set_ylabel('Hospitals')
# ax4[1, 0].set_ylabel('Hospitals')
#
# # Set x axis labels for bottom row
# ax4[1, 0].set_xlabel('Rooms')
# ax4[1, 1].set_xlabel('Rooms')
# ax4[1, 2].set_xlabel('Rooms')
#
# ax4[0, 0].set_title('Monday')
# ax4[0, 1].set_title('Tuesday')
# ax4[1, 0].set_title('Wednesday')
# ax4[1, 1].set_title('Thursday')
# ax4[1, 2].set_title('Friday')
#
# fig4.suptitle('Operating room cost per hour ($USD)')
#
# fig4.show()
#
# savename4 = 'Figures/K_hdr_heatmap.png'
# fig4.savefig(savename4, width=1980, height=1080)
