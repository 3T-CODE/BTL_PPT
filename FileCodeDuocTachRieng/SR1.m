% Copyright (C) 2024 Crypt
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {} {@var{retval} =} SR1 (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: Crypt <crypt@debian>
% Created: 2024-05-19

%SR1 is the command-line function:

% Khoi tao tinh toan bang SR1

function x_giaTri = SR1 (x_giaTri, f, e, delta, eta, r)
  global n
  global I
  x = sym ('x', [1, n]);

  % Khoi tao ma tran xap xi Hessian la ma tran don vi
  B = eye (n);

  % Tinh gradient ham so
  grad_f = gradient (f, x);
  grad_f_giaTri = double (subs (grad_f, x, x_giaTri));

  % Dat bien dem k = 0
  k = 0;

  % Khoi dong vong lap SR1
  while norm (grad_f_giaTri, "fro") > e
    % Tim s de q(s) = grad_f*s+0.5*s'*B*s dat gia tri nho nhat voi ||s|| <= delta
    s = TR_subproblem (f, delta, B, x_giaTri);

  % Tinh y , gia tri giam that , gia tri giam du doan
    % Tinh y
    grad_f_s_giaTri = double (subs (grad_f, x, x_giaTri + s));
    y = grad_f_s_giaTri - grad_f_giaTri;

    % Tinh gia tri giam that
    f_s_giaTri = double (subs (f, x, x_giaTri + s));
    f_giaTri = double (subs (f, x, x_giaTri));
    f_giam = f_giaTri - f_s_giaTri;

    % Tinh gia tri giam du doan
    f_giam_du_doan = -(grad_f_giaTri' * s + 0.5 * s' * B * s);

    % So sanh ty so f_giam/f_giam_du_doan voi eta
    if (f_giam / f_giam_du_doan) > eta
      x_giaTri = x_giaTri + s;
    else
      x_giaTri = x_giaTri;
    endif

    % So sanh ty so f_giam/f_giam_du_doan voi 0.75
    if (f_giam / f_giam_du_doan) > 0.75
      if norm (s, "fro") <= 0.8 * delta
        delta = delta;
      else
        delta = 2 * delta;
      endif
    elseif 0.1 <= (f_giam / f_giam_du_doan) && (f_giam / f_giam_du_doan) <= 0.75
      delta = delta;
    else
      delta = delta * 0.5;
    endif

    % So sanh quyet dinh cap nhat ma tran
    if abs (s' * (y - B * s)) >= r * norm (s, "fro") * norm (y - B * s, "fro")
      B = B + ((y - B * s) * (y - B * s)') / ((y - B * s)' * s);
    else
      B = B;
    endif

    % Tang bien k len 1 don vi
    k = k + 1;

    % Tinh gradient moi cua ham so
    grad_f_giaTri = double (subs (grad_f, x, x_giaTri));

  endwhile
endfunction

