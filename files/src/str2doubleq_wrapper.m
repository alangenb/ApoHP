function out = str2doubleq_wrapper(in)

out = real(str2doubleq(in));

% realized 2018-08-05 that str2doubleq() returns 0 instead of NaN
