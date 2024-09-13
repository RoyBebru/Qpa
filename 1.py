#!/usr/bin/env python

p = 3
s = 528259049275771151485822936813

# digits = []
# while True:
#     (s, r) = divmod(s, p)
#     digits.append(r)
#     if s == 0:
#         break

def loop_by_digits(s, p):
    while True:
        (s, r) = divmod(s, p)
        yield r
        if s == 0:
            break

def period_val_of(n: str):
  def period_of(n: str, _rolcheck = False):
    ln = len(n)
    if ln < 2:
      return n
    d = ln // 2
    while d > 0:
      nleft = n[:d]
      nright = n[d:d+d]
      if nleft == nright:

        if _rolcheck:
        # nleft=2200, tail=2200220022
        # Checking if nleft cover tail like nleft+nleft+nleft[:2]
          nright = n[d+d:] # tail
          dr = len(nright)
          if d >= dr:
          # nleft=2200, nright=22; nleft[:2] == nright?
            if nleft[:dr] != nright:
              return n
          else:
            for k in range(0, dr, d):
              if nleft != nright[k:k+d]:
                return n
            r = dr % d
            if r > 0 and nleft[-r:] != nright[k:]:
              return n

        return period_of(nleft, True)

      d -= 1
    return n # numbers without v like 10101101=(101)0

  w = period_of(n)
  lw = len(w)
  if lw == 0:
    return ('', n)
  ln = len(n)
  for k in range(0, ln, lw):
    v = n[k:k+lw]
    lv = len(v)
    if v != w:
      if lv > 0:
      # v is present
        # (2100)211 is the same as (0021)1
        for i in range(k, k+(lv if lv < lw else lw)):
          if w[0] != n[i]:
            break
          w = w[1:] + w[0]
        v = n[i:]
        break
    else:
      v = ''
  return (w, v)



digits = list(loop_by_digits(s, p))
digits.reverse()

print(digits)
s