
function accept=metropolis_log(d1,g1,d2,g2)

rat=exp(d2-d1)*exp(g2-g1);
if rat>1
   accept=1;
else
   r=rand;
   if r<rat
      accept=1;
   else
      accept=0;
   end
end
