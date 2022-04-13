#------------------------------------------------------------------------------------------------------------------------------------------#
##------------------------------------------------Routine for ELM sync and analysis-------------------------------------------------------##
#                                                     #-Guillermo Suarez Lopez-#                                                           #
#------------------------------------------------------------------------------------------------------------------------------------------#
#                                                       ##--Captain's log--##                                                              #
#------------------------------------------------------------------------------------------------------------------------------------------#
# Routines for ELM sync and analysis
#------------------------------------------------------------------------------------------------------------------------------------------#

##-- General Modules--##
import numpy as np
import os as os
import matplotlib.pyplot as plt
#import numba
import dd
from numpy.linalg import norm
from IPython import embed

##-- My Modules--##
#from gslpy.Data_Pkg import plot_arrange, plot_legend, plot_savefig, mov_aveg, close_to

#Packages to be imported
__all__ = ['_sync_grad','elm']

class elm():

    def __init__(self, shot=None, experiment='AUGD', ed=0, t1=0., t2=14., read_ELM=False):

        '''
        ----------------------------Description-----------------------------------
         Reads the ELM shotfile for ELM sync
        -------------------------------Input--------------------------------------
         read_ELM -> (Bool) Whether to read Mike's Dunne ELM shotfile.
                            Set to True to use this shotfile for ELM-sync.
        -------------------------------Output-------------------------------------
        ------------------------------Comments------------------------------------
        --------------------------------------------------------------------------
        '''

        self.shot       = shot
        self.experiment = experiment
        self.ed         = ed
	self.t1         = t1
	self.t2         = t2
	self.error_elm  = 0 

        elmc = {}
        if read_ELM:
            el               = dd.shotfile(diagnostic='ELM', pulseNumber=shot, experiment=experiment, edition=ed)
            t_begin          = el.getTimeBase('t_begELM', tBegin=t1, tEnd=t2)
            t_end            = el.getSignalCalibrated('t_endELM', tBegin=t1, tEnd=t2)[0]
            elmc['f_elm']    = el.getSignalCalibrated('f_ELM', tBegin=t1, tEnd=t2)[0]
            self.t_begin_elm = t_begin
            self.t_end_elm   = t_end
            self.elmc        = elmc
        else:
            self.error_elm = 1
            print('ACHTUNG: No ELM Shotfile, use Routine without ELM-related functions')

        return
#------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------#
 
    def sync_elm(self, timebase=None, signal=None, method=None, signalname='signal', 
                 time_preelm=None, time_elm=None, cycle=None, 
                 grad_thres_pos=None, grad_thres_neg=-600, safety_value=1., 
                 plot_flag=False):

        '''
        ----------------------------Description-----------------------------------
         Uses the "method" for ELM sync. Gives back tbegin and tend from ELMs
        -------------------------------Input--------------------------------------
         timebase -> (array) timebase of the signal to be synced
         signal   -> (array) Signal to be synced
         cycle    -> (list) Takes only the ELM cycle from 
                            {tend+cycle[0]*telm, tbeg-(1-cycle[1])*telm}
        -------------------------------Output-------------------------------------
	 timebase_elm  -> (array) timebase of only ELMs
        ------------------------------Comments------------------------------------
         The routine checks if it makes sense to ELM-sync in the first place.
            - If nu_s > 2*nu_elm, then the ELM can be well-resolved and it makes sense
            - If nu_s < 2*nu_elm, then ELMs are not resolved by the diagnostic
                                  and it makes no sense to ELM-sync.  

         The timebase is ELM-synced as: t_begin_i > t_j > t_end_i -> ELM
         with t_begin_i and t_end_i being found by method:
            - method = None      -> the ELM shotfile by M. Dunne.
            - method = Ipolsola  -> Gradient method on Ipolsola
            - method = Magnetics -> Using magnetic pick-up coils (NOT IMPLEMENTED)
        --------------------------------------------------------------------------
        '''

        #Checks if the sampling frequency of the signal is high enough (nu_s > 2*nu_elm)
        #that it makes sense to ELM-sync it in the first place
        if self.error_elm==0 and (len(timebase)/(timebase.max() - timebase.min())) < 2.*np.median(self.elmc['f_elm']):
            print('ACHTUNG: (nu_s = ' + str(len(timebase)/(timebase.max() - timebase.min())) +\
                  ' < 2.*median(nu_elm) - Doesnt make sense to ELMsync')
            return
        #Assume f_elm = 150 Hz
        else:
            if (len(timebase)/(timebase.max() - timebase.min())) < 2.*150:
                print('ACHTUNG: (nu_s = ' + str(len(timebase)/(timebase.max() - timebase.min())) +\
                      ' < 2.*median(150) - Doesnt make sense to ELMsync')
                return  

        #Use the ELM shotfile
        if method is None and self.error_elm == 1:
            print('None method cannot be chosen unless there is an "ELM" shotfile. Choose other method')
            return


        #Get {t_begin, t_end} from Ipolsola
        if method == 'Ipolsola':
            ma                 = dd.shotfile('MAC', self.shot)
            self.timebase_ipol = ma.getTimeBase('Ipolsola', dtype=np.float64, tBegin=self.t1, tEnd=self.t2)
            self.signal_ipol   = ma.getSignalCalibrated('Ipolsola', dtype=np.float64, tBegin=self.t1, tEnd=self.t2)

            #The IPOLSOLA signal is downsampled to generate 'stronger' gradients
            down_time = 20
            self.sync_grad(self.timebase_ipol[::down_time], self.signal_ipol[0][::down_time], grad_thres_pos=grad_thres_pos, grad_thres_neg=grad_thres_neg, 
                           time_preelm=time_preelm, time_elm=time_elm, safety_value=safety_value, signalname=signalname, plot_flag=plot_flag)

#        #NOT WORKING AT THE MOMENT
#        if method == 'Magnetics':
#            ma                   = dd.shotfile('MHA', self.shot)
#            self.timebase_B31_14 = ma.getTimeBase('B31-14', dtype=np.float64, tBegin=self.t1, tEnd=self.t2)
#            self.signal_B31_14     = ma.getSignalCalibrated('B31-14', dtype=np.float64, tBegin=self.t1, tEnd=self.t2)
#
#            self.sync_grad(self.timebase_B31_14, self.signal_B31_14, grad_thres_pos=grad_thres_pos, grad_thres_neg=grad_thres_neg, 
#                           time_preelm=time_preelm, time_elm=time_elm, safety_value=safety_value, signalname=signalname, plot_flag=plot_flag)


        #Modifies the beginning and ending times from the ELM cycle
        #Takes into account the specific duration between two contiguous ELMs
        if (cycle is not None) and (time_preelm is None) and (time_elm is None):
            t_elm   = np.abs(self.t_begin_elm[1:] - self.t_end_elm[:-1])   #ELM cycle length = end_i+1 - begin_i
            t_begin = np.hstack((self.t_begin_elm[0], self.t_begin_elm[1:] - (1-cycle[1])*t_elm))
            t_end   = np.hstack((self.t_end_elm[:-1] + cycle[0] * t_elm, self.t_end_elm[-1]))
        else:
            t_begin = self.t_begin_elm
            t_end   = self.t_end_elm

        #Computes the elm indexes of the signal
        if timebase is not None:
            timebase_matrix        = np.array([timebase,]*len(t_begin))
            t_begin_matrix         = np.array([t_begin,]*len(timebase)).T
            t_end_matrix           = np.array([t_end,]*len(timebase)).T
            bool_matrix            = np.multiply(timebase_matrix > t_begin_matrix, timebase_matrix < t_end_matrix)
            idxs_elm               = np.where(bool_matrix)[1]
            idxs_noelm             = np.setxor1d(np.arange(0, len(timebase)), idxs_elm)
            timebase_elmsync       = timebase[idxs_noelm]
            signal_elmsync         = signal[idxs_noelm]
            timebase_elm           = timebase[idxs_elm]
            signal_elm             = signal[idxs_elm]
            self.timebase_matrix   = timebase_matrix
            self.timebase_elm      = timebase_elm
            self.signal_elm        = signal_elm
            self.timebase_elmsync  = timebase_elmsync
            self.signal_elmsync    = signal_elmsync
            self.bool_matrix       = bool_matrix
            self.idxs_elm          = idxs_elm
            self.idxs_noelm        = idxs_noelm


        #If there is no indexes for ELM timetraces, inform
        if len(idxs_elm) == 0:
            print(' >>> NO ELMs were found <<< ')


        if plot_flag:
            plot_elm = plot_arrange('2d', 'normal', xlabel='Time [s]', ylabel='Signal [U]', title='', size='27',
                                    grid=True, nrows=2, ncols=1, sharex=True, sharey=True)
            plot_elm.ax[0][0].plot(timebase, signal, label=signalname)
            legend = plot_legend(plot_elm.ax[0][0])
            plot_elm.ax[0][0].plot(timebase_elmsync ,signal_elmsync,c='red',lw='2')
            plot_elm.ax[1][0].plot(timebase, signal,c='0.5')
            plot_elm.ax[1][0].plot(timebase_elmsync, signal_elmsync,lw='2',c='red',label='ELM Sync')
            legend = plot_legend(plot_elm.ax[1][0])
        
            for i in range(len(t_begin)):
                plot_elm.ax[0][0].plot((t_begin[i],t_begin[i]),(signal.min(),signal.max()),ls='--',c='black')
                plot_elm.ax[0][0].plot((t_end[i],t_end[i]),(signal.min(),signal.max()),ls='--',c='black')
#------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------#

    #Wrapper for numba
    def sync_grad(self, time, signal, grad_thres_pos=0.05, grad_thres_neg=-0.05, time_preelm=None, time_elm=None,  
                  safety_value=None, signalname=None, plot_flag=False):

       grad, elm_onset, idxs_clean, time_elmsync, signal_elmsync, t_begin_elm, t_end_elm = _sync_grad(time, signal,
		  grad_thres_pos=grad_thres_pos, grad_thres_neg=grad_thres_neg, time_preelm=time_preelm, time_elm=time_elm,  
                  safety_value=safety_value, signalname=signalname, plot_flag=plot_flag)
        
       self.grad           = grad
       self.elm_onset      = elm_onset
       self.idxs_clean     = idxs_clean
       self.time_elmsync   = time_elmsync
       self.signal_elmsync = signal_elmsync
       self.t_begin_elm    = t_begin_elm
       self.t_end_elm      = t_end_elm
#------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------#

#@numba.jit('Tuple((f8[:],i8[:],i8[:],f8[:],f8[:],f8[:],f8[:]))(f8[:],f8[:],f8,f8,f8,f8,f8,f8,b1)', nopython=False) 
#I have passed a f8 where the karg is a string (signalname). Idk how numba handles str
#@numba.jit
def _sync_grad(time, signal, grad_thres_pos=0.05, grad_thres_neg=-0.05, time_preelm=None, time_elm=None,  
               safety_value=None, signalname=None, method='gradient', plot_flag=False):

    '''
    ----------------------------Description-----------------------------------
     Uses different methods in order to get rid of ELMs
    -------------------------------Input--------------------------------------
    -------------------------------Output-------------------------------------
     t_begin_elm -> (Array) Beginning time of the ELMs
     t_end_elm   -> (Array) Ending time of the ELMs
    ------------------------------Comments------------------------------------
     method == 'gradient'
       - Gets rid of the ELMs by assuming a threshold on the signal gradient.
       - Recognizes ELM onsets by assuming the previous neighbour to be
         below the threshold and the next neighbour to be above it.
    --------------------------------------------------------------------------
    '''

#    embed()
    if method == 'gradient':

        fs = len(time)/(time[-1]-time[0])

        #Takes the signal gradient
        grad          = np.gradient(signal,1)
        idxs_array    = np.arange(0, len(grad))
        points_elm    = int(time_elm*fs)
        points_preelm = int(time_preelm*fs)

        print('   - ELM length = ' + str(time_preelm) + ', ' + str(time_elm) + ' s' + 
              ', which means ' + str(points_preelm+points_elm) + ' points ' + 'for fs=' + str(fs) + ' Hz')

        #Gets the ELM onset (characterized by grad > c*thres_up)
        #Get ELM onset whith the condition that the previous point is within grad range and the next one is not
        if grad_thres_pos is not None:
            grad_onset = np.where(grad > safety_value*grad_thres_pos)[0]
            grad_onset_prev = grad_onset - 1
            grad_onset_next = grad_onset + 1
            if grad_onset_next[-1] >= len(time):
                grad_onset_next[-1] = len(time) - 1
            prev = grad[grad_onset_prev] < safety_value*grad_thres_pos
            next = grad[grad_onset_next] > safety_value*grad_thres_pos

        else:
            grad_onset      = np.where(grad < safety_value*grad_thres_neg)[0]
            grad_onset_prev = grad_onset - 1
            grad_onset_next = grad_onset + 1
            if grad_onset_next[-1] >= len(time):
                grad_onset_next[-1] = len(time) - 1
            prev = grad[grad_onset_prev] > safety_value*grad_thres_neg
            next = grad[grad_onset_next] < safety_value*grad_thres_neg

        tot  = prev*next
        elm_onset = grad_onset[tot]


        #----------------------------------------------------------#
        # Now I add the additional condition that two ELM-onsets
        # should NOT be very close to each other.
        # I will use the arbitrary value of points_elm/4
        # Advantadges: Sampling frequency resilient
        # Disadvantadges: Arbitrary
        #----------------------------------------------------------#

        offset = elm_onset - points_elm/4.
        liers  = np.where(offset[1:] < elm_onset[:-1])[0]+1 #+1 because offset starts at [1:]
        elm_onset = np.delete(elm_onset, liers)
        #----------------------------------------------------------#


        #Beginning/end of the ELM
        elm_preelm    = elm_onset - points_preelm
        elm_endelm    = elm_onset + points_elm

        idxs_wrong = np.where(elm_endelm >= len(time))
        elm_endelm[idxs_wrong] = len(time) - 1


        #Remove any elm onset too close to the previous one, such that onset(i) + time_elm > onset(i+1)
        #PROBLEM: ELM very close together, on inter-ELM structures do not fulfill well the requirement.
        #I ended up with inter-ELM structures being removed and the next-ELM being accepted.

        #idx_norm      = idxs_array[grad_onset][1:] - idxs_array[grad_onset][0:-1]
        #elm_onset     = np.where(idx_norm > points_elm/4.)
        #elm_onset     = idxs_array[grad_onset[elm_onset]]
        #elm_onset     = np.hstack((idxs_array[grad_onset[0]], elm_onset))
        #test_onset    = elm_onset + points_elm
        #rest_test     = elm_onset[1:] - test_onset[:-1]
        #liers         = np.where(rest_test < 0)
        #elm_onset     = np.delete(elm_onset, liers)

        #Any point separated by less of telm from the previous gets discarded
        #rest_test    = idxs_array[grad_onset][1:] - idxs_array[grad_onset][0:-1]
        #liers        = np.where(rest_test < points_elm)[0]+1
        #print(points_elm)
        #print(liers)
        #elm_onset    = np.delete(idxs_array[grad_onset], liers)
        #print(idxs_array[grad_onset])
        #print(idxs_array[grad_onset][liers[1:]])
        #plt.figure()
        #plt.plot(rest_test)
        #print(elm_onset)


        #Gets the ELM and not ELM indexes
        idxs_matrix       = np.array([idxs_array,]*len(elm_onset))
        idxs_begin_matrix = np.array([elm_preelm,]*len(idxs_array)).T
        idxs_end_matrix   = np.array([elm_endelm,]*len(idxs_array)).T
        bool_matrix       = np.multiply(idxs_matrix > idxs_begin_matrix , idxs_matrix < idxs_end_matrix)
        idxs_noclean      = np.where(bool_matrix)[1]
        idxs_clean        = np.setxor1d(idxs_array, idxs_noclean)

        time_elmsync      = time[idxs_clean]
        signal_elmsync    = signal[idxs_clean]
        t_begin_elm       = time[elm_preelm]
        t_end_elm         = time[elm_endelm]

    
        if plot_flag:

            plt.figure()
            ax1 = plt.subplot(311)
            plt.title(signalname)
            plt.plot(time, grad, c='0.5')
            if grad_thres_pos is not None:
                plt.plot((time.min(), time.max()), (grad_thres_pos,grad_thres_pos), ls='--', c='black')
                plt.plot((time.min(), time.max()), (safety_value*grad_thres_pos,safety_value*grad_thres_pos), ls='--', c='red')
            else:
                plt.plot((time.min(), time.max()), (grad_thres_neg,grad_thres_neg), ls='--', c='black')
                plt.plot((time.min(), time.max()), (safety_value*grad_thres_neg,safety_value*grad_thres_neg), ls='--', c='red')
            plt.scatter(time[idxs_clean], grad[idxs_clean], c='red', s=150., label='Clean points')
            plt.scatter(time[elm_onset], grad[elm_onset], c='green', s=200., label='Onsets')
            plt.legend()

            ax2 = plt.subplot(312, sharex=ax1)
            plt.title(signalname)
            plt.plot(time, signal, c='0.5')
            plt.plot(time[idxs_clean], signal[idxs_clean], c='red')
            plt.scatter(t_begin_elm, signal[elm_preelm], c='magenta', s=150., label=r'$t_{being}^{elm}$')
            plt.scatter(t_end_elm, signal[elm_endelm], c='blue', s=150., label=r'$t_{end}^{elm}$')
            plt.legend()

            ax3 = plt.subplot(313, sharex=ax1)
            plt.scatter(time[idxs_clean], signal[idxs_clean], c='0.25')
            plt.xlim(time.min(),time.max())
            plt.ylim(signal[idxs_clean].min(), signal[idxs_clean].max())

    
    return grad, elm_onset, idxs_clean, time_elmsync, signal_elmsync, t_begin_elm, t_end_elm
#------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------#
