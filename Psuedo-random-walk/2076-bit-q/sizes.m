load "present_group.m";
r:=4;
w:=2^(32);
u:=Z!SquareRoot(w);
d:=Ceiling(Log(u,p-1));
t_bar:=2^Ceiling(Log(2,d*u));
t:=Z!(w/t_bar);
r_bar:=Z!(t/r);
w_bar:=Ceiling(p/w);

