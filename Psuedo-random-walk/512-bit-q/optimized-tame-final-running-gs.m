
load "present_group.m";
printf "\n q=%o",q;
load "target_bases.m";
printf "\n g1=%o",g1;
printf "\n g2=%o",g2;
printf "\n h=%o",h;


printf "\n dis_pt_length=%o",dis_pt_length;


N:=q;

tl:=0;
tb:=N;



load "table.m";

tame_seq:=[];
exp_tame_seq:=[];




pre_calculated_table_length:=#table;
printf "\n pre_calculated_table_length=%o",pre_calculated_table_length;

t_dis:=0; 


tame_initiations:=0;

sum_tame_walk_length:=0;



ii:=0;

load "fun.m";
d1:=dis_pt_length;tmp_itr_gs:=2^d1;
printf "\n tmp_itr_gs=%o",tmp_itr_gs;
pow_2_pre_calculated_table_length:= Z!Log(2,pre_calculated_table_length);

printf "\n time for %o tame points in GS is",length_for_tame_or_wild;
time while((#tame_seq lt length_for_tame_or_wild) ) do
	x1:=Random(tl,tb);
	x2:=Random(tl,tb);
	z:=Zp!((g1^x1)*(g2^x2));
	z_in_integers:=Z!z;
	distinguished_point:=ModByPowerOf2(z_in_integers,dis_pt_length) eq 0;	
	tame_initiations:=tame_initiations+1;
	length_of_tame_walk_now:=0;
	event_tame:=true;
	ii:=ii+1;
	while((not distinguished_point) and event_tame) do
		random_elm_index:=ModByPowerOf2(z_in_integers,pow_2_pre_calculated_table_length)+1;
		ii:=ii+1;
		random_elm_from_table:=table[random_elm_index];
		m:=random_elm_from_table[1];
		z:=Zp!prod_mod_p(z,m);	z_in_integers:=Z!z;
		x1:=x1+random_elm_from_table[2];
		x2:=x2+random_elm_from_table[3];			
		
		event_1:=(tl le x1) and (x1 le tb); 
		event_2:=(tl le x2) and (x2 le tb); 
		event_tame:=event_1 and event_2;
		length_of_tame_walk_now:=length_of_tame_walk_now+1;
		if(not event_tame) then
			sum_tame_walk_length:=sum_tame_walk_length+length_of_tame_walk_now;
			break;	
		else
			distinguished_point:=ModByPowerOf2(z_in_integers,dis_pt_length) eq 0;	
		end if;	
			
	end while;
	if(distinguished_point) then
		t_dis:=t_dis+1;
		Append(~tame_seq,z);
		Append(~exp_tame_seq,[x1,x2]);	
		//break;
	end if;
end while;

printf "\n points_walked_through=%o, points/(2^dis_pt_length)=%o",ii,RealField(5)!(ii/(2^dis_pt_length));

printf "\n points/(tame_length*tmp_itr_gs)=%o",RealField(5)!(ii/(tame_length*tmp_itr_gs));

printf "\n t_dis=%o, tame_initiations=%o, sum_tame_walk_length=%o",t_dis,tame_initiations,sum_tame_walk_length;
printf "\n number of distinct elements in tame seq=%o",#SequenceToSet(tame_seq);
