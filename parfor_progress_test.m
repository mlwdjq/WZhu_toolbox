N = 100;
parfor_progress(N);
for i=1:N
    pause(rand); % Replace with real code
    parfor_progress;
end
parfor_progress(0);