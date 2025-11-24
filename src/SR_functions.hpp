
#ifndef SR_FUNCTIONS_H
#define SR_FUNCTIONS_H


// --- Beverton–Holt: theta = (R0, h, SSB0)
template<class Type>
inline Type SR_BH(const vector<Type>& theta, const Type& SSB) {
  const Type R0   = theta(0);
  const Type h    = theta(1);
  const Type SSB0 = theta(2);
  const Type denom = (5.0 * h - 1.0);
  const Type alpha = (4.0 * h * R0) / denom;
  const Type beta  = SSB0 * (1.0 - h) / denom;
  return alpha * SSB / (beta + SSB);
}

// --- Ricker: theta = (alpha, beta)
template<class Type>
inline Type SR_Ricker(const vector<Type>& theta, const Type& SSB) {
  return theta(0) * SSB * exp(-theta(1) * SSB);
}

// --- Hockey-stick: theta = (alpha, beta)
template<class Type>
inline Type SR_Hockey(const vector<Type>& theta, const Type& SSB) {
  const Type alpha = theta(0);
  const Type beta  = theta(1);
  // AD-safe conditional
  return CppAD::CondExpLe(SSB, beta, alpha * SSB, alpha * beta);
}

// --- Dispatcher: recmodel 0=BH, 1=Ricker, 2=Hockey
template<class Type>
inline Type SR_eval(const int recmodel, const vector<Type>& theta, const Type& SSB) {
  switch (recmodel) {
    case 3: return SR_BH<Type>(theta, SSB);
    case 2: return SR_Ricker<Type>(theta, SSB);
    case 1: return SR_Hockey<Type>(theta, SSB);
    default: return SR_BH<Type>(theta, SSB);
  }
}

#endif
