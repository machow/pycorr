#import pieman_intersubj
#from scipy import io as sio
#from numpy.testing import assert_almost_equal
#M = sio.loadmat("janice/pieman_a1_3mm_intact_threesubs_roicorr.mat")
#j_sub = M['roitc']
#m_sub = np.array([tc / tc.std() for tc in roi_tcs]).T   #Our matrices are transposes
#m_sub.shape
#
#((j_sub - m_sub)**2).sum(axis=0) 	    #subs are aligned
#((j_sub[:,1] - m_sub[:,0])**2).sum(axis=0)  #ex of non-aligned
#
#m_ISC = intersubcorr(crosscor(m_sub.T))
#j_ISC = M['roicorr']
#np.assert_almost_equal(m_ISC, j_ISC)
