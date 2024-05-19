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
% @deftypefn {} {@var{retval} =} DFP (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: Crypt <crypt@debian>
% Created: 2024-05-19

%DFP is the command-line function:

% Khoi tao tinh toan bang DFP

function x_giaTri = DFP (x_giaTri, f, e)
  global n
  global I
  x = sym ('x', [1, n]);

  % Khoi tao ma tran xap xi Hessian dao la ma tran don vi
  H = eye (n);

  % Tinh gradient ham so
  grad_f = gradient (f, x);
  grad_f_giaTri = double (subs (grad_f, x, x_giaTri));

  % Vong lap xet gradient
  k = int8 (0);

  % Khoi dong vong lap DFP
  while norm (grad_f_giaTri, "fro") > e && k <= 10

  % Tim gia tri huong giam p
    p = -H * grad_f_giaTri;

    % Tim gia tri do dai buoc a thoa dieu kien Wolfe
    a = line_search (f, grad_f, x_giaTri, p);

    % Cap nhat x
    x_giaTri = x_giaTri + a * p;

    % Tinh s
    s = a * p;

    % Tinh y
    grad_f_giaTricu = grad_f_giaTri;
    grad_f_giaTri = double (subs (grad_f, x, x_giaTri));
    y = grad_f_giaTri - grad_f_giaTricu;

    % Cap nhat H
    H = H - (H * y * y' * H) / (y' * H * y) + (s * s') / (s' * y);

    % Cap nhat bien dem k
   k = k + 1;

    endwhile
endfunction

