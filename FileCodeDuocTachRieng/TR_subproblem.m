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
% @deftypefn {} {@var{retval} =} TR_subproblem (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: Crypt <crypt@debian>
% Created: 2024-05-19

%Khoi tao ham tim s thoa man bai toan phu
function p =  TR_subproblem( f , delta , B , x_giaTri)

   global n;
   x = sym('x', [1 n]);

   %Tinh gradient cua ham so
   grad_f = gradient(f,x);
   grad_f_giaTri = double( subs( grad_f , x ,  x_giaTri) );

   %Tinh p_s
   norm_grad_f_giaTri = norm(grad_f_giaTri , "fro");
   p_s = - (delta / norm_grad_f_giaTri) * grad_f_giaTri ;

   %Tinh gBg
   gBg = grad_f_giaTri' * B * grad_f_giaTri;

   %Xet gia tri cua gBg de tim tau
    if gBg <= 0
      tau = 1;

    elseif( norm_grad_f_giaTri.^3 / ( delta * grad_f_giaTri' * B * grad_f_giaTri ) < 1  )
        tau = norm_grad_f_giaTri.^3/( delta * grad_f_giaTri' * B * grad_f_giaTri )

      else
        tau = 1;

    end

    p = p_s * tau;

end

