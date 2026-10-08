function xy = xyzs2slepregionfixture(varargin)
    % Synthetic ocean-like region: an outer boundary, an island hole, and a
    % disconnected component. Accept GeoDomain's legacy name-value interface.
    if nargin == 1 && strcmp(varargin{1}, 'rotated')
        xy = false;
        return
    end

    xy = [20 -30; 70 -30; 65 35; 25 30; 20 -30; NaN NaN; ...
              35 -5; 35 10; 45 10; 45 -5; 35 -5; NaN NaN; ...
              90 0; 100 0; 100 15; 90 15; 90 0];
end
