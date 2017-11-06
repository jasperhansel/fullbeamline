!+
! Subroutine linear_fit (x, y, n_data, a, b, sig_a, sig_b)
!
! Subroutine to fit to y = A + B x
!-

subroutine linear_fit (x, y, n_data, a, b, sig_a, sig_b)

  use precision_def

  implicit none

  integer n, n_data

  real(rp) x(*), y(*), a, b, sig_a, sig_b
  real(rp) sum_x, sum_y, sum_xy, sum_x2, sum_y2, det
  real(rp) s2

  sum_x  = 0.0
  sum_y  = 0.0
  sum_xy = 0.0
  sum_x2 = 0.0
  sum_y2 = 0.0

  do n = 1, n_data
    sum_x  = sum_x  + x(n)
    sum_y  = sum_y  + y(n)
    sum_xy = sum_xy + x(n) * y(n)
    sum_x2 = sum_x2 + x(n)**2
    sum_y2 = sum_y2 + y(n)**2
  enddo

  det = n_data * sum_x2 - sum_x**2
  a = (sum_x2 * sum_y - sum_x * sum_xy) / det
  b = (n_data * sum_xy - sum_x * sum_y) / det

  s2 = (sum_y2 + a*a*n_data + b*b*sum_x2 +  &
            2*(a*b*sum_x - a*sum_y - b*sum_xy)) / n_data
  sig_a = sqrt(max(0.0_rp, s2 * sum_y2 / det))
  sig_b = sqrt(max(0.0_rp, s2 * n_data / det))

end subroutine

