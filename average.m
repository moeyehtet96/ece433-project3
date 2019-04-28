function avg = average(x,T,dt)
n_period = T/dt;
avg = 0;

for m = 1:n_period
    avg = avg + x(m)*dt/T;
end