import matplotlib.pyplot as plt
import matplotlib

import settings

def plot_counts_vs_time(selected_runs, roi_list):
    """
    Function to plot the 'Total Counts' value as a function of 'Start Time' for selected runs.

    Parameters:
        selected_runs (DataFrame): DataFrame containing the selected runs.
        roi_list (array): Array containing ROI definitions, ROI counts column header, ROI counts err column header
    """

    # Establish colorscheme
    cmap = matplotlib.colormaps['RdYlGn']

    for roi, roi_str, roi_err_str in roi_list:
        plt.figure(figsize=(10, 6))
        for i, B_value in enumerate(settings.B_select):
            runs_for_B = selected_runs[selected_runs['Field'] == B_value]
            duration = runs_for_B['Duration [Hours]']
            midpoint_time = runs_for_B['Start Time Elapsed [Hours]'] + duration/ 2

            #total_counts, total_counts_err = runs_for_B['Total Counts'], np.sqrt(runs_for_B['Total Counts'])
            counts, counts_err = runs_for_B[roi_str], runs_for_B[roi_err_str]

            color = cmap(i / len(settings.B_select))  # Sample color from colormap
            #plt.errorbar(midpoint_time, counts, xerr=half_duration, yerr=counts_err,
            #       capsize=2, ms=3, color=color, fmt=' ', label=f'|B| = {B_value}')
            plt.scatter(midpoint_time, counts,
                    s=30, color=color, label=f'|B| = {B_value}')

        plt.xlabel('Start Time [Hrs]')
        plt.ylabel(f'{roi_str} Total Counts')
        plt.title('Counts vs. Start Time')
        plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

        plt.tight_layout()
        plt.savefig(f'{settings.output_dir}{roi_str}_counts_plot.pdf', format='pdf')
        plt.show()

def plot_countrates_vs_time(selected_runs, roi_list):
    """
    Function to plot the 'Total Counts' value as a function of 'Start Time' for selected runs.

    Parameters:
        selected_runs (DataFrame): DataFrame containing the selected runs.
        roi_list (array): Array containing ROI definitions, ROI counts column header, ROI counts err column header
    """

    # Establish colorscheme
    cmap = matplotlib.colormaps['RdYlGn']

    for roi, roi_str, roi_err_str in roi_list:
        plt.figure(figsize=(10, 6))
        for i, B_value in enumerate(settings.B_select):
            runs_for_B = selected_runs[selected_runs['Field'] == B_value]
            duration = runs_for_B['Duration [Hours]']
            midpoint_time = runs_for_B['Start Time Elapsed [Hours]'] + duration/ 2

            countrates, countrates_err = runs_for_B[roi_str]/ duration, runs_for_B[roi_err_str]/ duration

            color = cmap(i / len(settings.B_select))  # Sample color from colormap
            plt.scatter(midpoint_time, countrates,
                    s=30, color=color, label=f'|B| = {B_value}')

        plt.xlabel('Start Time [Hrs]')
        plt.ylabel(roi_str + ' Count Rate [$hours^{-1}$]')
        plt.title('Count Rate vs. Start Time')
        plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

        plt.tight_layout()
        plt.savefig(f'{settings.output_dir}{roi_str}_countRate_plot.pdf', format='pdf')
        plt.show()

