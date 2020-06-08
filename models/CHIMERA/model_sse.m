% CHIMERA: Combining Mechanistic Models and Machine Learning for Personalized Chemotherapy and Surgery Sequencing in Breast Cancer
% Sum Squared Error, SSE
function sse = model_sse(alfa, sigma, M, y)
N = length(M);
estd = error_std(alfa, sigma, M, y);
sse = 0;
for id=1:N
    sse = sse + ((y(id) - M(id))/(estd(id) / sigma))^2;
end
end