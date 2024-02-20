q:=(1846389521368) + (11^600);//bit_q=2076
k:=2;
p:=k*q+1;

//Defining the group and the subgroup Gq 
Z:=Integers();
Zp := Integers(p);
Zp_star, iso_unit_group_to_Zp := UnitGroup(Zp);
Gq, iso_Gq_to_Zq := sub<Zp_star| ((p-1) div q)*Zp_star.1>;
length_for_tame_or_wild:=5;
tame_length:=length_for_tame_or_wild;
wild_length:=length_for_tame_or_wild;
dis_pt_length:=12;
dis_pt_length_tt:=6;
r1:=4;l1:=20;
Zq:=Integers(q);

