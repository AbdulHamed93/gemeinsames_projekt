% Definition der Wegpunkte
%3-Punkte Kurve
 waypoints = [0, 0;
              0, 5;
              6, 5;]

%4-Punkte Pfad
 % waypoints = [0, 0;
 %              0, 10;
 %              6, 10;
 %              6,0;]

 %5-Punkte Pfad
 waypoints = [3, 15;
              3,  5;
              10, 5;
              17, 5;
              17, 15;]

%Korridor-Punkte Pfad
% waypoints = [0, 0;
%              0, 1;
%              0, 2;
%              0, 3;
%              0, 4;
%              0, 5;
%              1, 5;
%              2, 5;
%              3, 5;
%              4, 5;
%              5, 5;
%              6, 5;
%              6, 4;
%              6, 3;
%              6, 2;
%              6, 1;
%              6,0;]

% Plotte die Kurven und äquidistante Punkte
step = 0.01; % Schrittweite für äquidistante Punkte
bezier_curve = computeBezierCurve(waypoints, step);
bspline_curve = computeBSplineCurve(waypoints, step);

figure;
plot(bezier_curve(:,1), bezier_curve(:,2), 'b-', 'LineWidth', 2); % Bezierkurve
hold on;
plot(bspline_curve(:,1), bspline_curve(:,2), 'r-', 'LineWidth', 2); % B-Spline-Kurve
plot(bezier_curve(:,1), bezier_curve(:,2), 'bx'); % Punkte auf Bezierkurve
plot(bspline_curve(:,1), bspline_curve(:,2), 'rx'); % Punkte auf B-Spline-Kurve
plot(waypoints(:,1), waypoints(:,2), 'blacko'); % Wegpunkte
xlabel('X');
ylabel('Y');
title('Kurven mit äquidistanten Punkten');
legend('Bezierkurve', 'B-Spline-Kurve','Zwischenpunkte Bezierkurve','Zwischenpunkte B-Spline-Kurve', 'Wegpunkte','Location', 'SouthEast');
grid on;

% Funktion zur Berechnung der Bezier-Kurve mit äquidistanten Zwischenpunkten
function bezier_curve = computeBezierCurve(waypoints, step)
    % Berechnung der Kontrollpunkte für die Bezier-Kurve
    control_points = zeros(size(waypoints, 1), 2);
    control_points(1,:) = waypoints(1,:);
    control_points(end,:) = waypoints(end,:);
    for i = 2:size(waypoints, 1)-1
        control_points(i,:) = (waypoints(i-1,:) + 2*waypoints(i,:) + waypoints(i+1,:)) / 4;
    end

    % Berechnung der Bezier-Kurve mit äquidistanten Zwischenpunkten
    t_bezier = 0:step:1;
    bezier_curve = zeros(length(t_bezier), 2);
    for i = 1:length(t_bezier)
        for j = 1:size(waypoints, 1)
            bezier_curve(i,:) = bezier_curve(i,:) + nchoosek(size(waypoints, 1)-1, j-1) * (1-t_bezier(i))^(size(waypoints, 1)-j) * t_bezier(i)^(j-1) * control_points(j,:);
        end
    end
end

% Funktion zur Berechnung der B-Spline-Kurve mit äquidistanten Zwischenpunkten
function bspline_curve = computeBSplineCurve(waypoints, step)
    % Berechnung der B-Spline-Kurve mit äquidistanten Zwischenpunkten
    degree = 2; % Grad des B-Splines
    pieces = size(waypoints, 1) - degree;
    knots = augknt(linspace(0, 1, pieces+1), degree+1);
    sp_x = spap2(knots, degree+1, linspace(0, 1, size(waypoints, 1)), waypoints(:,1));
    sp_y = spap2(knots, degree+1, linspace(0, 1, size(waypoints, 1)), waypoints(:,2));
    t_bspline = 0:step:1;
    bspline_curve = zeros(length(t_bspline), 2);
    for i = 1:length(t_bspline)
        bspline_curve(i,:) = [fnval(sp_x, t_bspline(i)), fnval(sp_y, t_bspline(i))];
    end
end