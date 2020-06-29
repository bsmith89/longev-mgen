import pymc3 as pm
import theano.tensor as tt
import theano
import scipy as sp
import numpy as np
from pymc3.distributions.continuous import PositiveContinuous
from pymc3.distributions.dist_math import bound
from pymc3.distributions.distribution import draw_values
import pymc3.distributions.special


def as_col(x):
    if isinstance(x, tt.TensorVariable):
        return x.dimshuffle(0, 'x')
    else:
        return np.asarray(x).reshape(-1, 1)

def as_row(x):
    if isinstance(x, tt.TensorVariable):
        return x.dimshuffle('x', 0)
    else:
        return np.asarray(x).reshape(1, -1)

def exp(x):
    if isinstance(x, tt.TensorVariable):
        return tt.exp(x)
    else:
        return np.exp(x)

def log(x):
    if isinstance(x, tt.TensorVariable):
        return tt.log(x)
    else:
        return np.log(x)

def _ifthen(pred, iftrue, iffalse):
    assert (len(pred) == len(iftrue)) and (len(pred) == len(iffalse))
    out = np.empty_like(pred)
    for i, (p, t, f) in enumerate(zip(pred, iftrue, iffalse)):
        if p:
            out[i] = t
        else:
            out[i] = f
    return out

def switch(x, iftrue, iffalse):
    if isinstance(x, tt.TensorVariable):
        return tt.switch(x, iftrue, iffalse)
    else:
        return _ifthen(x, iftrue, iffalse)

def lt(x, y):
    if isinstance(x, tt.TensorVariable) or isinstance(y, tt.TensorVariable):
        return tt.lt(x, y)
    else:
        return np.less(x, y)

def dot(x, y):
    if isinstance(x, tt.TensorVariable) or isinstance(y, tt.TensorVariable):
        return tt.dot(x, y)
    else:
        return np.dot(x, y)

def gammaln(x):
    if isinstance(x, (tt.Constant, tt.Variable)):
        return pymc3.distributions.special.gammaln(x)
    else:
        return sp.special.gammaln(x)

def sum(xx, *args, **kwargs):
    if isinstance(xx, (tt.TensorVariable, tt.TensorType,
                       tt.TensorConstant, np.ndarray)):
        return xx.sum(*args, **kwargs)
    else:
        if args or kwargs:
            raise ValueError("args or kwargs was defined, "
                             "but xx is not an array type.")
        else:
            return sum(xx)

def dm_logp(x, alpha, n=None):
    if not n:
        n = sum(x, axis=-1)
    sum_alpha = sum(alpha, axis=-1)

    const = (gammaln(n + 1) + gammaln(sum_alpha)) - gammaln(n + sum_alpha)
    series = gammaln(x + alpha) - (gammaln(x + 1) + gammaln(alpha))
    result = const + sum(series, axis=-1)
    return result


class DirichletMultinomial(pm.Discrete):
    """Dirichlet Multinomial Distribution

    Parameters
    ----------
    alpha : tensor
        dirichlet parameter

    """
    def __init__(self, alpha, n, *args, **kwargs):
        super(DirichletMultinomial, self).__init__(*args, **kwargs)
        self.alpha = alpha
        self.n = n

    def logp(self, x):
        alpha = self.alpha
        return dm_logp(x, alpha)

    def random(self, point=None, size=None, repeat=None):
        raise NotImplementedError
        # alpha, n = draw_values([self.alpha, self.n], point=point)
        # return rv_dm(alpha=alpha, n=n, size=size)


def rv_d(alpha, size=None):
    """Generate samples from the Dirichlet distribution."""
    alpha = np.asarray(alpha)
    k = alpha.shape[-1]
    if size is not None:
        assert alpha.ndim == 1, "Can not specify size parameter with a 2-d alpha."
        s = size
        alpha = np.outer(np.ones(s), alpha)
    elif alpha.ndim == 2:
        s = alpha.shape[0]
    else:
        s = 1
        alpha = np.outer(np.ones(s), alpha)

    out = np.empty((s, k))
    for i in range(s):
        out[i, :] = sp.random.dirichlet(alpha[i,:])
    if size is None and s == 1:
        return out[0,:]
    else:
        return out


def rv_dm(alpha, n, size=None):
    """Generate samples from the Dirichlet-Multinomial distribution."""
    alpha = np.asarray(alpha)
    k = alpha.shape[-1]
    n = np.asarray(n)
    if size is not None:
        assert alpha.ndim == 1 and n.ndim == 0, "Can not specify size parameter with a 2-d alpha."
        s = size
        alpha = np.outer(np.ones(s), alpha)
        n = np.outer(np.ones(s), n)
    elif alpha.ndim == 2:
        s = alpha.shape[0]
        if n.ndim == 1:
            assert alpha.shape[0] == n.shape[0], "Length of *n* array does not match length of *alpha* array"
        else:
            n = np.outer(np.ones(s), n)
    elif n.ndim == 1:
        s = alpha.shape[0]
        alpha = np.outer(np.ones(s), alpha)
    else:
        s = 1
        n = np.outer(np.ones(s), n)
        alpha = np.outer(np.ones(s), alpha)

    out = np.empty((s, k))
    for i in range(s):
        out[i, :] = sp.random.multinomial(n[i], sp.random.dirichlet(alpha[i,:]))
    if size is None and s == 1:
        return out[0,:]
    else:
        return out


def norm_simplex(p):
    """Sum-to-zero transformation."""
    return (p.T / sum(p, axis=-1)).T

def backwards_interval(x, a=0, b=1):
    """Backwards interval transformation.

    Maps from the whole real line to the interval [0, 1].

    """
    return (b - a) * np.exp(x) / (1 + np.exp(x)) + a

def backwards_stickbreaking(p):
    """Backwards stickbreaking transformation.

    Maps from a vector of length (k - 1)
    to a k-simplex coordinate.

    """
    assert p.ndim == 1
    Km1 = p.shape[0]
    k = np.arange(Km1)[(slice(None), ) + (None, ) * (p.ndim - 1)]
    eq_share = -np.log(Km1 - k)
    z = 1 / (1 + np.exp(-(p + eq_share)))
    p1 = np.concatenate([z, np.ones(p[:1].shape)])
    pu = np.concatenate([np.ones(p[:1].shape), 1 - z])
    S = np.cumprod(pu, 0)
    x = S * p1
    return x

def ccmodel(beta, x):
    """Community composition model."""
    return norm_simplex(exp(dot(x, log(beta))))

def along_env(comm, pert, env):
    """Generate community composition vector with a proportional perturbation.

    Takes 1-d vectors for pert/env instead of 2-d.

    See :along_envs:.

    TODO: Implement using ccmodel().

    """
    pert = switch(lt(pert.ndim, 2), as_row(pert), pert)
    env = switch(lt(env.ndim, 2), as_col(env), env)
    return along_envs(comm, pert, env)

def along_envs(comm, pert, env):
    """Generate community composition vectors with proportional perturbations.

    Specifically, the perturbation is applied as the product of
    the base community and the perturbation vectors raised to the power of
    the environmental factor.

    comm : tensor (1-d)
        An order-$k$ row-vector
        (treated as though it is on
        the $(k-1)$ unit simplex)

    pert : tensor (2-d)
        An order $p \\times k$ matrix
        where each row is on the $k-1$ unit
        simplex (sums to 1)
        ($p$ environmental parameters
        $k$ species)

    env : tensor (2-d)
        An order $n \\times p$ matrix
        ($n$ samples by $p$ environmental
        parameters)

    $$
    \\frac{\\exp(X \\dot \\log(\\Delta) + \log(\\pi))}
         {\\sum{\\exp(X \\dot \\log(\\Delta) + \\log(\\pi))}}
    $$

    Where $X$ is the environmental vector, $\\Delta$ is
    the environmental effect matrix, and $\\pi$ is the
    basal community (under no environmental impact).

    TODO: Implement using ccmodel().

    """
    return norm_simplex(exp(dot(env, log(pert)) + log(comm)))

class OffsetLognormal(PositiveContinuous):
    def __init__(self, mu=0, tau=1, offset=0, *args, **kwargs):
        """Lognormal distribution with some missed detections."""
        super(OffsetLognormal, self).__init__(*args, **kwargs)
        self.mu = mu
        self.tau = tau
        self.offset = offset
        self.mean = tt.exp(self.mu) - self.offset

    def logp(self, value):
        mu = self.mu
        tau = self.tau
        offset = self.offset
        v = value + offset
        log_pdf = (-0.5 * tau * (tt.log(v) - mu)**2
                   + 0.5 * tt.log(tau / (2. * np.pi))
                   - tt.log(v))
        return bound(log_pdf, tau > 0, offset > 0)
