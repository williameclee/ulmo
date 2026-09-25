function epochAveragingTest
    % Synthetic data checks independent of external product files.
    dates = datetime(2020, 1, 1) + days(0:59)';
    mesh = reshape(1:60, 1, 1, []);
    epochs = [datetime(2020, 1, 1), datetime(2020, 2, 1); ...
                  datetime(2020, 2, 1), datetime(2020, 3, 1); ...
                  datetime(2020, 4, 1), datetime(2020, 5, 1)];
    [values, labels] = averageepochs(dates, mesh, epochs);
    assert(isequal(squeeze(values(1, 1, 1:2)), [16; 46]));
    assert(isnan(values(1, 1, 3)));
    assert(isequal(labels, mean(epochs, 2)));
    mesh(1, 1, 1) = NaN;
    values = averageepochs(datenum(dates), mesh, datenum(epochs));
    assert(values(1, 1, 1) == 16.5);
    dates = datetime(2020, 1, 1) + days(0:99)';
    mesh = cat(1, reshape(1:100, 1, 1, []), reshape(101:200, 1, 1, []));
    [values, labels] = interptemporal(dates, mesh, days(10), 'linear', BeQuiet = true);
    [~, group] = min(abs(dates(:)' - labels(:)), [], 1);

    for k = 1:numel(labels)
        assert(isequal(values(:, :, k), mean(mesh(:, :, group == k), 3)));
    end

    reference = [datetime(2017, 6, 3), datetime(2017, 6, 28, 18, 0, 0); ...
                     datetime(2018, 6, 5, 6, 0, 0), datetime(2018, 6, 29)];
    [filled, synthetic] = fillgraceepochs(reference, [datetime(2017, 6, 1), datetime(2018, 7, 1)]);
    gap = filled(synthetic, :);
    assert(gap(1, 1) == reference(1, 2));
    assert(gap(end, 2) == reference(2, 1));
    assert(isequal(gap(1:end - 1, 2), gap(2:end, 1)));
    assert(isequal(filled(~synthetic, :), reference));
    [subset, ~] = fillgraceepochs(reference, datetime(2017, [9 11], 16));
    labels = mean(filled, 2);
    keep = labels >= datetime(2017, 9, 16) & labels <= datetime(2017, 11, 16);
    assert(isequal(subset, filled(keep, :)));
    mustBeTimeStep('grace');
    disp('epochAveragingTest passed');
end
