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
% @deftypefn {} {@var{retval} =} line_search (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: Crypt <crypt@debian>
% Created: 2024-05-19

%Khoi tao ham tim do dai buoc a
function a =  line_search( f , grad_f , x_giaTri , p   )
  global n;
  x = sym('x', [1 n]);
  c1 = 0.0001;
  c2 = 0.9;
  a = 1 ;

  f_a_giaTri = double(subs(f , x , x_giaTri + a*p));
  f_giaTri = double(subs(f , x , x_giaTri));
  grad_f_giaTri = double(subs(grad_f , x , x_giaTri));
  grad_f_a_giaTri = double(subs(grad_f , x , x_giaTri + a * p));

 %Xet a theo dieu kien Wolfe ||  abs( ((-1) * grad_f_a_giaTri' * p) ) > abs( ((-1) * c2 * grad_f_giaTri' * p) )
  while  f_a_giaTri > f_giaTri + c1 * a * grad_f_giaTri' * p ||  abs( ((-1) * grad_f_a_giaTri' * p) ) > abs( ((-1) * c2 * grad_f_giaTri' * p) )


                a = 0.5*a;
                f_a_giaTri = double(subs(f , x , x_giaTri + a*p));
                f_giaTri = double(subs(f , x , x_giaTri));
                grad_f_giaTri = double(subs(grad_f , x , x_giaTri));
                grad_f_a_giaTri = double(subs(grad_f , x , x_giaTri + a * p));


  end

end



