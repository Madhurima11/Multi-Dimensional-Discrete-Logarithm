q:=670390396497129854978701249910292306373968291029619668886178072186088201503\
6773488400937149083451713845015929093243025426876941405973284973216824506305919;//bit_q=512
k:=2;
p:=k*q+1;

//Defining the group and the subgroup Gq 
Z:=Integers();
Zp := Integers(p);
Zp_star, iso_unit_group_to_Zp := UnitGroup(Zp);
Gq, iso_Gq_to_Zq := sub<Zp_star| ((p-1) div q)*Zp_star.1>;
length_for_tame_or_wild:=100;
tame_length:=length_for_tame_or_wild;
wild_length:=length_for_tame_or_wild;
dis_pt_length:=10;
dis_pt_length_tt:=5;
r1:=4;l1:=20;
Zq:=Integers(q);

