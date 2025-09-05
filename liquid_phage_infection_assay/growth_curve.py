import os
import pandas as pd
import numpy as np
from datetime import datetime, time, timedelta
import matplotlib.pyplot as plt
import scipy.stats
import seaborn as sns
from sklearn import metrics

#%%
def fit_growth_splines(time_series,logOD_series,spline_size=15,r_sq_threshold=0.98,trust_OD_range=[0.05,0.8]):
    # return best linear regression fits for spline_size data points, as well as saturation OD

    # calculate saturation
    saturation_logOD = np.max(logOD_series)
    saturation_OD = np.exp(saturation_logOD)

    max_found = False
    # break into segments
    for t1 in range(len(time_series)):
        t2 = t1+spline_size
        if t2>len(time_series):
            #end of dataframe
            continue
        t = time_series[t1:t2]
        y = logOD_series[t1:t2]

        # check for nans, untrustworthy logOD values
        trust_logOD_range = np.log(trust_OD_range)
        if not(y.between(trust_logOD_range[0],trust_logOD_range[1]).any()):
            # contains nans
            continue

        # calculate fit
        fit = scipy.stats.linregress(t, y)

        if (fit.rvalue**2>=r_sq_threshold):
            if not(max_found):
                best_fit = fit
                best_t = t
                max_found = True

            if fit.slope>best_fit.slope:
                best_fit = fit
                best_t = t
                max_found = True
    if max_found:
        return best_fit, best_t, saturation_OD
    else:
        return None, None, saturation_OD

def calc_area_under_curve(time_series):
    auc = np.trapz(time_series, axis=0)
    return auc

def plot_growth_curve(time_series,logOD_series,plot_title,best_fit=None,best_t=None,saturation_OD=None,save_plots=True):

    # Plot Paramteters
    plot_ylim = [-5,1]
    plot_xlim = [0,np.ceil(time_series.iat[-1])]
    #plot_xlim = [0,24]
    plot_xlabel = 'Time [h]'
    plot_ylabel ='log(OD)'

    fig,ax = plt.subplots()
    ax.plot(time_series,logOD_series,color='black', linestyle='None', marker='o', markersize=2)
    ax.set(xlabel=plot_xlabel, ylabel=plot_ylabel,title=plot_title)
    ax.set_ylim(plot_ylim)
    ax.set_xlim(plot_xlim)

    if saturation_OD:
        ax.plot(plot_xlim,[np.log(saturation_OD),np.log(saturation_OD)],'b',linewidth=1)

    if best_fit:
        ax.plot(best_t, best_fit.intercept + best_fit.slope*best_t, 'r', linewidth=1)
        ax.text(best_t.iloc[0]+1,best_fit.intercept + best_fit.slope*best_t.iloc[0]-0.5,'GR = ' + str(round(best_fit.slope,3)) + ' 1/h')
        ax.text(best_t.iloc[0]+1,best_fit.intercept + best_fit.slope*best_t.iloc[0]-1,'R^2 = ' + str(round(best_fit.rvalue,3)))

    fig.savefig('Figures/' + plot_title + '.png')
    return ax

def process_plate(file_name,plate_number,curvefit_parameters):
    plate_name = 'Plate ' + str(plate_number) +' - Raw Data'
    data = pd.read_excel(file_name,plate_name,skiprows=[0],dtype={'Time': datetime})

    #more complicated version of what was originally here because the first code couldn't read past 24 hours
    raw_time_series = data['Time']
    elapsed_time = []
    og_time = datetime(year = 1899, month = 12, day = 31, hour = 0, minute = 0, second = 0)
    fix_time_dict = {}

    for t in raw_time_series:
        if type(t) == time:
            string = str(t)
            list = string.split(sep = ":")
            date_time = datetime(year = 1899, month = 12, day = 31, hour = int(list[0]), minute = int(list[1]), second = int(list[2]))
        else:
            date_time = t
        time_delta = date_time - og_time
        elapsed_hours = ((time_delta.days)*24 + (time_delta.seconds)/3600)
        elapsed_time.append(elapsed_hours)
        fix_time_dict[t] = elapsed_hours
    time_series = pd.Series(elapsed_time)

    #correct data
    ODdata = data.iloc[:,1:]
    ODdata = ODdata.apply(lambda x: np.log(x-curvefit_parameters['blank']))
    ODdata = ODdata.rolling(curvefit_parameters['roll_size']).mean()
    
    results = pd.DataFrame(index=range(len(ODdata.columns)), columns=['Well','Growth Rate', 'Time to Max Growth', 'Saturation OD', 'Rsq'])

    for n in range(len(ODdata.columns)):
        plot_title = 'Plate ' + str(plate_number) + ' Well ' + ODdata.columns[n]
        plot_series = ODdata.iloc[:,n]
        best_fit,best_t,saturation_OD = fit_growth_splines(time_series,plot_series,curvefit_parameters['spline_size'],r_sq_threshold=0.98)
        
        if curvefit_parameters['save_plots']:
            gc_ax = plot_growth_curve(time_series,plot_series,plot_title,best_fit,best_t,saturation_OD)

        if best_fit:
            growth_rate =  best_fit.slope
            time_to_max = best_t.iloc[0]
            rsq = best_fit.rvalue
        else:
            growth_rate = 0
            time_to_max = 0
            rsq = 0
        
        
        results.loc[[n],:] = pd.DataFrame(index = [n],data={'Well':plot_title, 'Growth Rate':growth_rate, 'Time to Max Growth':time_to_max, 'Saturation OD':saturation_OD, 'Rsq':rsq})

        ODdata_hours = data.set_index('Time').rename(index = fix_time_dict)
    return results, ODdata_hours

#%%
if __name__ == "__main__":

    #user parameters
    file_name = 'platesAandBcombination.xlsx'
    output_name = 'platesAandBcombination'
    out_directory = 'Data/'
    plate_numbers = [1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    auc_time_bound = 72

    blank = 0.09
    roll_size = 15  #rolling average over n values
    spline_size = 2*roll_size+1 #calculate growth rate over segments of n values
    curvefit_parameters = {'blank':blank, 'roll_size':roll_size, 'spline_size':spline_size, 'save_plots':True}

    # matplotlib parameters
    plt.close('all')
    plt.rcParams['figure.max_open_warning'] = 0

    # make folder for figures
    if not(os.path.exists('Figures')):
        os.mkdir('Figures')
        os.mkdir('Data')
    allresults = pd.DataFrame(columns=['Well','Growth Rate', 'Time to Max Growth', 'Saturation OD', 'Rsq', 'AUC'])
    all_ODs = {}
    for plate_number in plate_numbers:

        print('Processing Plate ' +str(plate_number) + ' of ' + str(len(plate_numbers)))
        results_list = process_plate(file_name,plate_number,curvefit_parameters) 
        
        # all_ODs[plate_number] = results_list[1]
        
        tmp_OD = results_list[1].rolling(curvefit_parameters['roll_size']).mean().dropna()
        all_ODs[plate_number] = tmp_OD
        
        tmp_results = results_list[0]
        tmp_results['AUC'] = np.trapz(tmp_OD.iloc[tmp_OD.index <= auc_time_bound,:], axis=0)
        
        allresults = pd.concat([allresults,tmp_results],ignore_index=True)
        
    allresults.to_csv(out_directory+output_name + '_results.csv')
    
    #save od timeseries from each plate
    for plate_index in all_ODs:
        pd.DataFrame(all_ODs[plate_index]).to_csv(out_directory+"all_ODs_plate"+str(plate_index)+".csv")



# %% 



# %% 





# %% 





