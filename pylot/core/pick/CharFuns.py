# -*- coding: utf-8 -*-
"""
Created Oct/Nov 2014

Implementation of the Characteristic Functions (CF) published and described in:

Kueperkoch, L., Meier, T., Lee, J., Friederich, W., & EGELADOS Working Group, 2010:
Automated determination of P-phase arrival times at regional and local distances
using higher order statistics, Geophys. J. Int., 181, 1159-1170

Kueperkoch, L., Meier, T., Bruestle, A., Lee, J., Friederich, W., & EGELADOS
Working Group, 2012: Automated determination of S-phase arrival times using
autoregressive prediction: application ot local and regional distances, Geophys. J. Int.,
188, 687-702.

:author: MAGS2 EP3 working group
"""
import numpy as np
from obspy.core import Stream

class CharacteristicFunction(object):
    '''
    SuperClass for different types of characteristic functions.
    '''
    def __init__(self, data, cut, t2=None, order=None, t1=None, fnoise=0.001):
        '''
        Initialize data type object with information from the original
        Seismogram.

        :param: data
        :type: `~obspy.core.stream.Stream`

        :param: cut
        :type: tuple

        :param: t2
        :type: float

        :param: order
        :type: int

        :param: t1
        :type: float (optional, only for AR)

        :param: fnoise
        :type: float (optional, only for AR)
        '''

        assert isinstance(data, Stream), "%s is not a stream object" % str(data)

        self.orig_data = data
        self.dt = self.orig_data[0].stats.delta
        self.setCut(cut)
        self.setTime1(t1)
        self.setTime2(t2)
        self.setOrder(order)
        self.setFnoise(fnoise)
        self.setARdetStep(t2)
        self.calcCF(self.getDataArray())
        self.arpara = np.array([])
        self.xpred = np.array([])

    def __str__(self):
        return '''\n\t{name} object:\n
        Cut:\t\t{cut}\n
        t1:\t{t1}\n
        t2:\t{t2}\n
        Order:\t\t{order}\n
        Fnoise:\t{fnoise}\n
        ARdetStep:\t{ardetstep}\n
        '''.format(name=type(self).__name__,
                   cut=self.getCut(),
                   t1=self.getTime1(),
                   t2=self.getTime2(),
                   order=self.getOrder(),
                   fnoise=self.getFnoise(),
                   ardetstep=self.getARdetStep[0]())

    def getCut(self):
        return self.cut

    def setCut(self, cut):
        self.cut = cut

    def getTime1(self):
        return self.t1

    def setTime1(self, t1):
        self.t1 = t1

    def getTime2(self):
        return self.t2

    def setTime2(self, t2):
        self.t2 = t2

    def getARdetStep(self):
        return self.ARdetStep

    def setARdetStep(self, t1):
        if t1:
           self.ARdetStep = []
           self.ARdetStep.append(t1 / 4)
           self.ARdetStep.append(int(np.ceil(self.getTime2() / self.getIncrement()) / 4))

    def getOrder(self):
        return self.order

    def setOrder(self, order):
        self.order = order

    def getIncrement(self):
        """
        :rtype : int
        """
        return self.dt

    def getTimeArray(self):
        if self.getTime1():
           incr = self.getARdetStep()
           self.TimeArray = np.arange(0, len(self.getCF()) * incr[0], incr[0]) + self.getCut()[0] \
                           + self.getTime1() + self.getTime2()  
        else:
           incr = self.getIncrement()
           self.TimeArray = np.arange(0, len(self.getCF()) * incr, incr) + self.getCut()[0]  
        return self.TimeArray

    def getFnoise(self):
        return self.fnoise

    def setFnoise(self, fnoise):
        self.fnoise = fnoise

    def getCF(self):
        return self.cf

    def getDataArray(self, cut=None):
        '''
        If cut times are given, time series is cut from cut[0] (start time)
        till cut[1] (stop time) in order to calculate CF for certain part
        only where you expect the signal!
        input: cut (tuple) ()
        cutting window
        '''
        if cut is not None:
            if self.cut[0] == 0:
                start = 0
            else:
                start = self.cut[0] / self.dt
                stop = self.cut[1] / self.dt
            if len(self.orig_data) == 1:
                zz = self.orig_data.copy()
                z1 = zz[0].copy()
                zz[0].data = z1.data[start:stop]
                data = zz
                return data
            elif len(self.orig_data) == 2:
                hh = self.orig_data.copy()
                h1 = hh[0].copy()
                h2 = hh[1].copy()
                hh[0].data = h1.data[start:stop]
                hh[1].data = h2.data[start:stop]
                data = hh 
                return data
            elif len(self.orig_data) == 3:
                hh = self.orig_data.copy()
                h1 = hh[0].copy()
                h2 = hh[1].copy()
                h3 = hh[2].copy()
                hh[0].data = h1.data[start:stop]
                hh[1].data = h2.data[start:stop]
                hh[2].data = h3.data[start:stop]
                data = hh 
                return data
        else:
            data = self.orig_data.copy()
            return data
       
    def calcCF(self, data=None):
        self.cf = data


class AICcf(CharacteristicFunction):
    '''
    Function to calculate the Akaike Information Criterion (AIC) after
    Maeda (1985).
    :param: data, time series (whether seismogram or CF)
    :type: tuple

    Output: AIC function
    '''

    def calcCF(self, data):

        print 'Calculating AIC ...'
        x = self.getDataArray()
        xnp = x[0].data
        datlen = len(xnp)
        k = np.arange(1, datlen)
        cf = np.zeros(datlen)
        cumsumcf = np.cumsum(np.power(xnp, 2))
        i = np.where(cumsumcf == 0)
        cumsumcf[i] = np.finfo(np.float64).eps
        cf[k] = ((k - 1) * np.log(cumsumcf[k] / k) + (datlen - k + 1) *
                 np.log((cumsumcf[datlen - 1] -
                        cumsumcf[k - 1]) / (datlen - k + 1)))
        cf[0] = cf[1]
        inf = np.isinf(cf)
        ff = np.where(inf == True)
        if len(ff) >= 1:
            cf[ff] = 0

        self.cf = cf - np.mean(cf)


class HOScf(CharacteristicFunction):
    '''
    Function to calculate skewness (statistics of order 3) or kurtosis
    (statistics of order 4), using one long moving window, as published
    in Kueperkoch et al. (2010).
    '''

    def calcCF(self, data):

        x = self.getDataArray(self.getCut())
        xnp =x[0].data
        if self.getOrder() == 3:  # this is skewness
            print 'Calculating skewness ...'
            y = np.power(xnp, 3)
            y1 = np.power(xnp, 2)
        elif self.getOrder() == 4:  # this is kurtosis
            print 'Calculating kurtosis ...'
            y = np.power(xnp, 4)
            y1 = np.power(xnp, 2)

        #Initialisation
        #t2: long term moving window
        ilta = round(self.getTime2() / self.getIncrement())
        lta = y[0]
        lta1 = y1[0]

        #moving windows
        LTA = np.zeros(len(xnp))
        for j in range(0, len(xnp)):
            if j < 4:
                LTA[j] = 0
            elif j <= ilta:
                lta = (y[j] + lta * (j-1)) / j
                lta1 = (y1[j] + lta1 * (j-1)) / j
            else:
                lta = (y[j] - y[j - ilta]) / ilta + lta
                lta1 = (y1[j] - y1[j - ilta]) / ilta + lta1
            #define LTA
            if self.getOrder() == 3:
                LTA[j] = lta / np.power(lta1, 1.5)
            elif self.getOrder() == 4:
                LTA[j] = lta / np.power(lta1, 2)
        
        self.cf = LTA


class ARZcf(CharacteristicFunction):

    def calcCF(self, data):

        print 'Calculating AR-prediction error from single trace ...'
        x = self.getDataArray(self.getCut())
        xnp = x[0].data
        #some parameters needed
        #add noise to time series
        xnoise = xnp + np.random.normal(0.0, 1.0, len(xnp)) * self.getFnoise() * max(abs(xnp))
        tend = len(xnp)
        #Time1: length of AR-determination window [sec]
        #Time2: length of AR-prediction window [sec]
        ldet = int(round(self.getTime1() / self.getIncrement()))    #length of AR-determination window [samples]
        lpred = int(np.ceil(self.getTime2() / self.getIncrement())) #length of AR-prediction window [samples]

        cf = []
        loopstep = self.getARdetStep()
        for i in range(ldet + self.getOrder() - 1, tend - lpred + 1, loopstep[1]):
            #determination of AR coefficients
            self.arDetZ(xnoise, self.getOrder(), i-ldet, i)
            #AR prediction of waveform using calculated AR coefficients
            self.arPredZ(xnp, self.arpara, i + 1, lpred)
            #prediction error = CF
            err = np.sqrt(np.sum(np.power(self.xpred[i:i + lpred] - xnp[i:i + lpred], 2)) / lpred)
            cf.append(err)
        
        #convert list to numpy array
        cf = np.asarray(cf)
        self.cf = cf


    def arDetZ(self, data, order, rind, ldet):
        '''
        Function to calculate AR parameters arpara after Thomas Meier (CAU), published
        in  Kueperkoch et al. (2012). This function solves SLE using the Moore-
        Penrose inverse, i.e. the least-squares approach.
        :param: data, time series to calculate AR parameters from
        :type: array

        :param: order, order of AR process
        :type: int

        :param: rind, first running summation index
        :type: int

        :param: ldet, length of AR-determination window (=end of summation index)
        :type: int

        Output: AR parameters arpara
        '''

        #recursive calculation of data vector (right part of eq. 6.5 in Kueperkoch et al. (2012)
        rhs = np.zeros(self.getOrder())
        for k in range(0, self.getOrder()):
            for i in range(rind, ldet):
                rhs[k] = rhs[k] + data[i] * data[i - k]

        #recursive calculation of data array (second sum at left part of eq. 6.5 in Kueperkoch et al. 2012)
        A = np.zeros((2,2))
        for k in range(1, self.getOrder() + 1):
            for j in range(1, k + 1):
                for i in range(rind, ldet):
                    ki = k - 1
                    ji = j - 1
                    A[ki,ji] = A[ki,ji] + data[i - ji] * data[i - ki]

                A[ji,ki] = A[ki,ji]

        #apply Moore-Penrose inverse for SVD yielding the AR-parameters
        self.arpara = np.dot(np.linalg.pinv(A), rhs)

    def arPredZ(self, data, arpara, rind, lpred):
        '''
        Function to predict waveform, assuming an autoregressive process of order
        p (=size(arpara)), with AR parameters arpara calculated in arDet. After 
        Thomas Meier (CAU), published in Kueperkoch et al. (2012).
        :param: data, time series to be predicted
        :type: array

        :param: arpara, AR parameters
        :type: float

        :param: rind, first running summation index
        :type: int

        :param: lpred, length of prediction window (=end of summation index)
        :type: int

        Output: predicted waveform z
        '''
        #be sure of the summation indeces
        if rind < len(arpara) + 1:
            rind = len(arpara) + 1
        if rind > len(data) - lpred + 1:
            rind = len(data) - lpred + 1
        if lpred < 1:
            lpred = 1
        if lpred > len(data) - 1:
            lpred = len(data) - 1

        z = np.append(data[0:rind], np.zeros(lpred))
        for i in range(rind, rind + lpred):
            for j in range(1, len(arpara) + 1):
                ji = j - 1
                z[i] = z[i] + arpara[ji] * z[i - ji]

        self.xpred = z


class ARHcf(CharacteristicFunction):

    def calcCF(self, data):

        print 'Calculating AR-prediction error from both horizontal traces ...'
        
        xnp = self.getDataArray(self.getCut())
        
        #some parameters needed
        #add noise to time series
        xenoise = xnp[0].data + np.random.normal(0.0, 1.0, len(xnp[0].data)) * self.getFnoise() * max(abs(xnp[0].data))
        xnnoise = xnp[1].data + np.random.normal(0.0, 1.0, len(xnp[1].data)) * self.getFnoise() * max(abs(xnp[1].data))
        Xnoise = np.array( [xenoise.tolist(), xnnoise.tolist()] )
        tend = len(xnp[0].data)
        #Time1: length of AR-determination window [sec]
        #Time2: length of AR-prediction window [sec]
        ldet = int(round(self.getTime1() / self.getIncrement()))    #length of AR-determination window [samples]
        lpred = int(np.ceil(self.getTime2() / self.getIncrement())) #length of AR-prediction window [samples]
              
        cf = []
        loopstep = self.getARdetStep()
        for i in range(ldet + self.getOrder() - 3, tend - lpred + 1, loopstep[1]):
            self.arDetH(Xnoise, self.getOrder(), i-ldet, i)
            #AR prediction of waveform using calculated AR coefficients
            self.arPredH(xnp, self.arpara, i + 1, lpred)
            #prediction error = CF
            err = np.sqrt(np.sum(np.power(self.xpred[0][i:i + lpred] - xnp[0][i:i + lpred], 2) \
            + np.power(self.xpred[1][i:i + lpred] - xnp[1][i:i + lpred], 2)) / (2 * lpred))
            cf.append(err)

        #convert list to numpy array
        cf = np.asarray(cf)
        self.cf = cf

    def arDetH(self, data, order, rind, ldet):
        '''
        Function to calculate AR parameters arpara after Thomas Meier (CAU), published
        in  Kueperkoch et al. (2012). This function solves SLE using the Moore-
        Penrose inverse, i.e. the least-squares approach. "data" is a structured array.
        AR parameters are calculated based on both horizontal components in order
        to account for polarization.
        :param: data, horizontal component seismograms to calculate AR parameters from
        :type: structured array

        :param: order, order of AR process
        :type: int

        :param: rind, first running summation index
        :type: int

        :param: ldet, length of AR-determination window (=end of summation index)
        :type: int

        Output: AR parameters arpara
        '''

        #recursive calculation of data vector (right part of eq. 6.5 in Kueperkoch et al. (2012)
        rhs = np.zeros(self.getOrder())
        for k in range(0, self.getOrder()):
            for i in range(rind, ldet):
                rhs[k] = rhs[k] + data[0,i] * data[0,i - k] + data[1,i] * data[1,i - k]

        #recursive calculation of data array (second sum at left part of eq. 6.5 in Kueperkoch et al. 2012)
        A = np.zeros((4,4))
        for k in range(1, self.getOrder() + 1):
            for j in range(1, k + 1):
                for i in range(rind, ldet):
                    ki = k - 1
                    ji = j - 1
                    A[ki,ji] = A[ki,ji] + data[0,i - ji] * data[0,i - ki] + data[1,i - ji] *data[1,i - ki]

                A[ji,ki] = A[ki,ji]

        #apply Moore-Penrose inverse for SVD yielding the AR-parameters
        self.arpara = np.dot(np.linalg.pinv(A), rhs)

    def arPredH(self, data, arpara, rind, lpred):
        '''
        Function to predict waveform, assuming an autoregressive process of order
        p (=size(arpara)), with AR parameters arpara calculated in arDet. After 
        Thomas Meier (CAU), published in Kueperkoch et al. (2012).
        :param: data, horizontal component seismograms to be predicted
        :type: structured array

        :param: arpara, AR parameters
        :type: float

        :param: rind, first running summation index
        :type: int

        :param: lpred, length of prediction window (=end of summation index)
        :type: int

        Output: predicted waveform z
        :type: structured array
        '''
        #be sure of the summation indeces
        if rind < len(arpara) + 1:
            rind = len(arpara) + 1
        if rind > len(data[0]) - lpred + 1:
            rind = len(data[0]) - lpred + 1
        if lpred < 1:
            lpred = 1
        if lpred > len(data[0]) - 1:
            lpred = len(data[0]) - 1

        z1 = np.append(data[0][0:rind], np.zeros(lpred))
        z2 = np.append(data[1][0:rind], np.zeros(lpred))
        for i in range(rind, rind + lpred):
            for j in range(1, len(arpara) + 1):
                ji = j - 1
                z1[i] = z1[i] + arpara[ji] * z1[i - ji]
                z2[i] = z2[i] + arpara[ji] * z2[i - ji]

        z = np.array( [z1.tolist(), z2.tolist()] )
        self.xpred = z

class AR3Ccf(CharacteristicFunction):

    def calcCF(self, data):

        print 'Calculating AR-prediction error from all 3 components ...'
        
        xnp = self.getDataArray(self.getCut())
        
        #some parameters needed
        #add noise to time series
        xenoise = xnp[0].data + np.random.normal(0.0, 1.0, len(xnp[0].data)) * self.getFnoise() * max(abs(xnp[0].data))
        xnnoise = xnp[1].data + np.random.normal(0.0, 1.0, len(xnp[1].data)) * self.getFnoise() * max(abs(xnp[1].data))
        xznoise = xnp[2].data + np.random.normal(0.0, 1.0, len(xnp[2].data)) * self.getFnoise() * max(abs(xnp[2].data))
        Xnoise = np.array( [xenoise.tolist(), xnnoise.tolist(), xznoise.tolist()] )
        tend = len(xnp[0].data)
        #Time1: length of AR-determination window [sec]
        #Time2: length of AR-prediction window [sec]
        ldet = int(round(self.getTime1() / self.getIncrement()))    #length of AR-determination window [samples]
        lpred = int(np.ceil(self.getTime2() / self.getIncrement())) #length of AR-prediction window [samples]
              
        cf = []
        loopstep = self.getARdetStep()
        for i in range(ldet + self.getOrder() - 3, tend - lpred + 1, loopstep[1]):
            self.arDet3C(Xnoise, self.getOrder(), i-ldet, i)
            #AR prediction of waveform using calculated AR coefficients
            self.arPred3C(xnp, self.arpara, i + 1, lpred)
            #prediction error = CF
            err = np.sqrt(np.sum(np.power(self.xpred[0][i:i + lpred] - xnp[0][i:i + lpred], 2) \
            + np.power(self.xpred[1][i:i + lpred] - xnp[1][i:i + lpred], 2) \
            + np.power(self.xpred[2][i:i + lpred] - xnp[2][i:i + lpred], 2)) / (3 * lpred))
            cf.append(err)

        #convert list to numpy array
        cf = np.asarray(cf)
        self.cf = cf

    def arDet3C(self, data, order, rind, ldet):
        '''
        Function to calculate AR parameters arpara after Thomas Meier (CAU), published
        in  Kueperkoch et al. (2012). This function solves SLE using the Moore-
        Penrose inverse, i.e. the least-squares approach. "data" is a structured array.
        AR parameters are calculated based on both horizontal components and vertical 
        componant.
        :param: data, horizontal component seismograms to calculate AR parameters from
        :type: structured array

        :param: order, order of AR process
        :type: int

        :param: rind, first running summation index
        :type: int

        :param: ldet, length of AR-determination window (=end of summation index)
        :type: int

        Output: AR parameters arpara
        '''

        #recursive calculation of data vector (right part of eq. 6.5 in Kueperkoch et al. (2012)
        rhs = np.zeros(self.getOrder())
        for k in range(0, self.getOrder()):
            for i in range(rind, ldet):
                rhs[k] = rhs[k] + data[0,i] * data[0,i - k] + data[1,i] * data[1,i - k] \
                + data[2,i] * data[2,i - k]

        #recursive calculation of data array (second sum at left part of eq. 6.5 in Kueperkoch et al. 2012)
        A = np.zeros((4,4))
        for k in range(1, self.getOrder() + 1):
            for j in range(1, k + 1):
                for i in range(rind, ldet):
                    ki = k - 1
                    ji = j - 1
                    A[ki,ji] = A[ki,ji] + data[0,i - ji] * data[0,i - ki] + data[1,i - ji] *data[1,i - ki] \
                    + data[2,i - ji] *data[2,i - ki]

                A[ji,ki] = A[ki,ji]

        #apply Moore-Penrose inverse for SVD yielding the AR-parameters
        self.arpara = np.dot(np.linalg.pinv(A), rhs)

    def arPred3C(self, data, arpara, rind, lpred):
        '''
        Function to predict waveform, assuming an autoregressive process of order
        p (=size(arpara)), with AR parameters arpara calculated in arDet3C. After 
        Thomas Meier (CAU), published in Kueperkoch et al. (2012).
        :param: data, horizontal and vertical component seismograms to be predicted
        :type: structured array

        :param: arpara, AR parameters
        :type: float

        :param: rind, first running summation index
        :type: int

        :param: lpred, length of prediction window (=end of summation index)
        :type: int

        Output: predicted waveform z
        :type: structured array
        '''
        #be sure of the summation indeces
        if rind < len(arpara) + 1:
            rind = len(arpara) + 1
        if rind > len(data[0]) - lpred + 1:
            rind = len(data[0]) - lpred + 1
        if lpred < 1:
            lpred = 1
        if lpred > len(data[0]) - 1:
            lpred = len(data[0]) - 1

        z1 = np.append(data[0][0:rind], np.zeros(lpred))
        z2 = np.append(data[1][0:rind], np.zeros(lpred))
        z3 = np.append(data[2][0:rind], np.zeros(lpred))
        for i in range(rind, rind + lpred):
            for j in range(1, len(arpara) + 1):
                ji = j - 1
                z1[i] = z1[i] + arpara[ji] * z1[i - ji]
                z2[i] = z2[i] + arpara[ji] * z2[i - ji]
                z3[i] = z3[i] + arpara[ji] * z3[i - ji]

        z = np.array( [z1.tolist(), z2.tolist(), z3.tolist()] )
        self.xpred = z
