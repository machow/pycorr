import numpy as np
from funcs_correlate import standardize, corsubs, lagcor
from numpy.testing import assert_almost_equal
import os

class test_standardize:
    """
    Test funcs_correlate.standardize
    """
    def setup(self):
        self.S1 = np.array([1.,2,3])

    def test_returns_copy_by_default(self):              #returns copy by default
        assert self.S1 is not standardize(self.S1)

    def test_returns_inplace_arg(self):                  #can work inplace
        tmp_S1  = self.S1.copy()
        assert tmp_S1 is standardize(tmp_S1, inplace=True)

    def test_inplace_equals_copy(self):                  #equal results for copy and inplace
        S1_inplace = standardize(self.S1, inplace=True)
        assert_almost_equal(standardize(self.S1), S1_inplace)

    def test_inplace_int(self):                          #int arrays never operate inplace
        tmp_S1 = self.S1.astype(int)
        assert tmp_S1 is not standardize(tmp_S1, inplace=True)
        assert_almost_equal(standardize(tmp_S1), standardize(tmp_S1, inplace=True))

    def test_axis_arg(self):                             #axis arg returns arr with expected shape
        tmp_M = np.array([np.ones(3), np.zeros(3)])
        tmp_stand = standardize(tmp_M, axis=0)
        assert_almost_equal(tmp_stand.shape, tmp_M.shape)

    def test_no_demean(self):
        print self.S1
        print standardize(self.S1, demean=False)
        assert not np.allclose(standardize(self.S1, demean=False).mean(), 0)

    def test_python_list(self):                         #can take python array
        el = list(self.S1)
        a = [el, el]
        M = np.array([self.S1, self.S1])
        assert_almost_equal(standardize(a), standardize(M))


class test_corsubs:
    """
    Test funcs_correlate.corsubs
    TODO: find better 3-dimensional array to test (anscombe only to 3 decimal prec
    """
    @classmethod
    def setup_class(self): 
        basedir = os.path.split(os.path.abspath(__file__))[0]
        fname_anscombe = os.path.join(basedir, 'data/anscombe.csv')
        #a, b, solution, decimal
        self.ts = {}
        anscomb = np.genfromtxt(fname_anscombe, delimiter=',', names=True)   #columns are x1, .., x4, y1, .., y4
        for ii in range(1,5):
            anscombe_xy = [anscomb['x%s'%ii], anscomb['y%s'%ii]]
            stand_xy = map(lambda s: standardize(s), anscombe_xy)
            sol = [.8165, 3]
            self.ts['anscombe%s'%ii] = anscombe_xy + sol
            self.ts['anscombe_stand%s'%ii] = stand_xy + sol

        a = np.arange(11)
        self.ts['linear'] = [a, a*2 + 10, 1]
        #make 3 dim array
        test = ['anscombe%s'%ii for ii in range(1,5)] + ['linear']
        #fix for order issues with dict
        self.a_3dim = np.array([v[0] for k, v in self.ts.iteritems() if k in test])
        self.b_3dim = np.array([v[1] for k, v in self.ts.iteritems() if k in test])
        self.corrs_3dim = np.array([v[2] for k, v in self.ts.iteritems() if k in test])

    def test_cor_expected_unstandardized(self):
        for k, v in self.ts.iteritems():
            yield self.cor_expected, k

    def test_cor_expected_standardized(self):
        for k, v in self.ts.iteritems():
            if 'stand' in k: 
                yield self.cor_expected, k, False

    def cor_expected(self, ts_key, standardize=False):
        x,y,sol = self.ts[ts_key][0:3]
        decimal = self.ts[ts_key][3] if len(self.ts[ts_key]) == 4 else 7
        corr = corsubs(x,y, standardized=standardize)
        assert_almost_equal(corr, sol, decimal=decimal)

    def test_lagcor_equal_subcor(self):
        for k, v in self.ts.iteritems():
            yield self.lagcor_equal_subcor, k

    def test_lagcor_equal_subcor_standardized(self):
        for k, v in self.ts.iteritems():
            if 'stan' in k:
                yield self.lagcor_equal_subcor, k, True

    def lagcor_equal_subcor(self, key, standardized=False):
        x, y = self.ts[key][0:2]
        assert_almost_equal(lagcor(x,y,0), corsubs(x,y))


    def test_cor_3d(self):
        assert_almost_equal(corsubs(self.a_3dim, self.b_3dim), self.corrs_3dim, decimal=3)

    def test_axis_arg(self):
        assert_almost_equal(corsubs(self.a_3dim.T, self.b_3dim.T, axis=0), self.corrs_3dim.T, decimal=3)

    def test_shape_mismatch(self):
        try: 
            corsubs(np.arange(10), np.arange(11))
            raise BaseException('returned correlation for mismatched dims')
        except ValueError: pass


    def test_crosscor_equal_subcor(self):
        #TODO implement between group crosscor
        pass


class test_lagcor:
    def setup(self):
        pass

    def test_lagcor_equal_shift(self):
        pass

