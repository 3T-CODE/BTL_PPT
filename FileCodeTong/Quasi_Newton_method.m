%3 Phuong phap tuong tu Newton
%Tac gia : Nhom 22 , Lop PPT L10

pkg load symbolic;

%Thiet lap lai gia tri
clear all;
close all;

%Nhap gia tri ban dau cua x , ham so f(x) , va gia tri dung sai hoi tu e
x_giaTri = input("Nhap gia tri x ban dau : x = ");
global n;
n = int16(length(x_giaTri));

sym f;
x = sym('x', [1 n]);
global I ;
I = eye(n);

f = input("Nhap ham so : f(x) = ");
e = double(input("Nhap gia tri dung sai hoi tu : e = "));

%Kiem tra gia tri cua e
while e <= 0
  e = double(input("Nhap gia tri dung sai hoi tu : e = "));
end





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





%Khoi tao tinh toan bang BFGS
function x_giaTri = BFGS( x_giaTri , f , e)

 global n;
 global I;
 x = sym('x', [1 n]);

 %Khoi tao ma tran xap xi Hessian dao la ma tran don vi
 H = eye(n);

 %Tinh gradient ham so
 grad_f = gradient(f,x);
 grad_f_giaTri = double(subs(grad_f , x , x_giaTri ));

 %Vong lap xet gradient
 k = int8(0);

 %Khoi dong vong lap BFGS
 while norm(grad_f_giaTri , "fro") > e && k <= 10

 %Tim gia tri huong giam p
        p = - H*grad_f_giaTri;

  %Tim gia tri do dai buoc a thoa dieu kien Wolfe
        a =line_search(f , grad_f , x_giaTri , p );

   %Cap nhat x
        x_giaTri = x_giaTri + a*p;

    %Tinh s
        s = a*p;

     %Tinh y
        grad_f_giaTricu = grad_f_giaTri;
        grad_f_giaTri = double(subs(grad_f , x , x_giaTri ));
        y = grad_f_giaTri - grad_f_giaTricu ;

    %Cap nhat H
        H = ( I - (s * y')/(y'*s) ) * H * ( I - (y*s')/(y'*s) ) + (s*s')/(y'*s) ;
        k=k+1;


 end

end





%Khoi tao tinh toan bang DFP
function x_giaTri = DFP( x_giaTri , f , e)

 global n;
 global I;
 x = sym('x', [1 n]);

 %Khoi tao ma tran xap xi Hessian dao la ma tran don vi
 H = eye(n);

 %Tinh gradient ham so
 grad_f = gradient(f,x);
 grad_f_giaTri = double(subs(grad_f , x , x_giaTri ));

 %Vong lap xet gradient
 k = int8(0);

 %Khoi dong vong lap DFP
 while norm(grad_f_giaTri , "fro") > e && k <= 10

 %Tim gia tri huong giam p
        p = - H*grad_f_giaTri;

  %Tim gia tri do dai buoc a thoa dieu kien Wolfe
        a =line_search(f , grad_f , x_giaTri , p );

   %Cap nhat x
        x_giaTri = x_giaTri + a*p;

    %Tinh s
        s = a*p;

     %Tinh y
        grad_f_giaTricu = grad_f_giaTri;
        grad_f_giaTri = double(subs(grad_f , x , x_giaTri ));
        y = grad_f_giaTri - grad_f_giaTricu ;

    %Cap nhat H
        H = H - (H * y * y' * H) / (y' * H * y) + (s * s') / (s' * y);
        k=k+1;


 end

end





%Khoi tao tinh toan bang SR1
function x_giaTri = SR1(x_giaTri , f , e  , delta , eta , r)
  global n
  global I
   x = sym('x', [1 n]);

   %Khoi tao ma tran xap xi Hessian la ma tran don vi
   B = eye(n);

   %Tinh gradient ham so
   grad_f = gradient( f , x );
   grad_f_giaTri = double(subs( grad_f , x , x_giaTri ));

   %Dat bien dem k = 0
   k=0 ;

   %Khoi dong vong lap SR1

   while norm(grad_f_giaTri ,"fro") > e



   %Tim s de q(s) = grad_f*s+0.5*s'*B*s dat gia tri nho nhat voi ||s|| <= delta
   s = TR_subproblem( f , delta , B , x_giaTri );

   %Tinh y , gia tri giam that , gia tri giam du doan
   %%Tinh y
   grad_f_s_giaTri = double( subs( grad_f , x , x_giaTri + s ));
   y = grad_f_s_giaTri  - grad_f_giaTri ;

   %%Tinh gia tri giam that
   f_s_giaTri = double( subs( f , x , x_giaTri + s ) );
   f_giaTri = double( subs( f , x , x_giaTri ) );
   f_giam = f_giaTri - f_s_giaTri ;

   %%Tinh gia tri giam du doan
   f_giam_du_doan = -(grad_f_giaTri' * s + 0.5* s' * B * s);

   %So sanh ty so f_giam/f_giam_du_doan voi eta
   if (f_giam / f_giam_du_doan) > eta
     x_giaTri = x_giaTri + s;

   else
     x_giaTri = x_giaTri ;

   end

   %So sanh ty so f_giam/f_giam_du_doan voi 0.75
   if (f_giam / f_giam_du_doan) > 0.75
        if norm(s , "fro") <= 0.8 * delta
                delta = delta;

         else
                delta = 2*delta;

        end

      elseif  0.1 <= (f_giam / f_giam_du_doan) && (f_giam / f_giam_du_doan) <= 0.75
        delta = delta ;

      else
        delta = delta*0.5 ;

      end


  %So sanh quyet dinh cap nhat ma tran
  if abs( s' * (y - B*s) ) >= r * norm(s , "fro") * norm( y - B*s,"fro")
        B = B +  ( ( y - B*s )* ( y - B*s )' ) / ( ( y - B*s )' * s );

      else
        B = B ;


 end

 %Tang bien k len 1 don vi
   k = k +1 ;

  %Tinh gradient moi cua ham so
   grad_f_giaTri = double(subs( grad_f , x , x_giaTri ));

 end

end




%Chon cach tinh
fprintf("\n Ban muon tim gi ? \n");
fprintf(" \n 1/Tim cuc tieu gan dung - (Nhap 1)");
fprintf(" \n 2/Tim cuc dai gan dung - (Nhap 2)");
fprintf(" \n 3/Tim nghiem gan dung - (Nhap 3) ");
fprintf(" \n 4/Huy tinh toan - (Nhap so khac) \n");
Lua_chon_cach_tinh = input("Nhap : ") ;

%Xet lua chon da nhap
if Lua_chon_cach_tinh == 1
  %Chuyen doi ham so phu hop
  f = f ;

  %Thong bao lua chon da nhap
  fprintf("\nBan da chon tim cuc tieu \n");

elseif Lua_chon_cach_tinh == 2
  %Chuyen doi ham so phu hop
  f = -f ;

  %Thong bao lua chon da nhap
  fprintf("\nBan da chon tim cuc dai \n");

elseif Lua_chon_cach_tinh == 3
  %Chuyen doi ham so phu hop
  f = f.^2;

  %Thong bao lua chon da nhap
  fprintf("\nBan da chon tim nghiem \n");

else

end

%Lua chon phuong phap tinh
fprintf("\n Hay lua chon phuong phap tinh ");
fprintf("\n 1/ BFGS - (Nhap 1)");
fprintf("\n 2/ DFP  - (Nhap 2)");
fprintf("\n 3/ SR1 - (Nhap 3) ");
fprintf("\n 4/ Huy tinh toan  - (Nhap so khac) \n");
Lua_chon_PP_tinh = input("Nhap : ");

%Xet lua chon da nhap va khoi dong phuong phap tinh da chon
if Lua_chon_PP_tinh == 1
  %Thong bao PP da lua chon
  fprintf("\n Ban da chon phuong phap BFGS \n ");

  %Khoi dong tinh toan
  x_giaTri = BFGS( x_giaTri , f , e);

  %In gia tri ra man hinh
  fprintf("\n Dap an : \n");
  fprintf("\n X = \n");
  disp(x_giaTri);

elseif Lua_chon_PP_tinh == 2
  %Thong bao PP da lua chon
  fprintf("\n Ban da chon phuong phap DFP \n ");

  %Khoi dong tinh toan
  x_giaTri = DFP( x_giaTri , f , e);

  %In gia tri ra man hinh
  fprintf("\n Dap an : \n");
  fprintf("\n X = \n");
  disp(x_giaTri);

elseif Lua_chon_PP_tinh == 3
   %Thong bao PP da chon
   fprintf("\n Ban da chon phuong phap SR1 \n");

   %Nhap cac du lieu con thieu
   fprintf("\n Hay nhap cac gia tri can thiet sau day : \n");
   delta = input("\n Nhap ban kinh hoi tu (Delta > 0):  Delta = ");
   eta = input("\n Nhap toan tu eta thuoc (0 ; 10^-3) : eta = ");
   r = input("\n Nhap toan tu r thuoc (0;1) : r = ");

   %Kiem tra du lieu da dat dieu kien hay chua
   while 1
        if delta > 0
          if 0 < eta && eta < 1e-3
            if 0 < r && r < 1
              break;
            end
          end
        end

       %Tien hanh nhap lai du lieu
    fprintf("\n Hay nhap cac gia tri can thiet sau day : \n");
   delta = input("\n Nhap ban kinh hoi tu (Delta > 0):  Delta = ");
   eta = input("\n Nhap toan tu eta thuoc (0 ; 10^-3) : eta = ");
   r = input("\n Nhap toan tu r thuoc (0;1) : r = ");

   end

   %Tim gia tri theo SR1
   x_giaTri = SR1(x_giaTri , f , e  , delta , eta , r);

   %In gia tri da tim duoc ra man hinh
   fprintf("\n Dap an : \n");
   fprintf("\nX =  \n");
   disp(x_giaTri);

else
   %Thong bao PP da nhap
   fprintf("\n Ban da huy gia tri tinh toan \n");


end
%%%%Ket thuc %%%%










