function metrics = validate_nqpca(classicalLambda,qpcaLambda)
%VALIDATE_NQPCA Compare N-QPCA/QPE estimates with classical reference.

classicalLambda = classicalLambda(:);
qpcaLambda = qpcaLambda(:);

K=min(numel(classicalLambda),numel(qpcaLambda));

absError=abs(classicalLambda(1:K)-qpcaLambda(1:K));
relError=absError./max(classicalLambda(1:K),eps);

metrics.absoluteError=absError;
metrics.relativeError=relError;
metrics.meanAbsoluteError=mean(absError);
metrics.maxAbsoluteError=max(absError);
metrics.topKClassicalExplainedPercent=100*sum(classicalLambda(1:K));
metrics.topKQPCAPhasePercent=100*sum(qpcaLambda(1:K));

end
