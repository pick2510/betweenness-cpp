//
// Created by std on 10.05.19.
//

#ifndef DEM_UTILS_ACCUMULATOR_H
#define DEM_UTILS_ACCUMULATOR_H

template <class T> class Accumulator {
private:
  T acc_val{0};

public:
  T getAccVal() const;

private:
  T max{0};
  T min{0};
  int n{0};

public:
  void operator()(T val);
  inline void update_min_max(T val);
  T get_mean() { return acc_val / n; };
};

template <class T> void Accumulator<T>::operator()(T val)
{
  acc_val += val;
  n++;
  update_min_max(val);
}

template <class T> inline void Accumulator<T>::update_min_max(T val)
{
  if (n == 1) {
    max = val;
    min = val;
  }
  else {
    val > max ? max = val : val;
    val < min ? min = val : val;
  }
}
template <class T> T Accumulator<T>::getAccVal() const { return acc_val; }

#endif // DEM_UTILS_ACCUMULATOR_H
