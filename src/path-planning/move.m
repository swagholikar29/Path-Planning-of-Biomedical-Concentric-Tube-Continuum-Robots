% TODO: direction needs to be multiple R3 unit vectors
function q_new = move(q_near, q_rand, delta_q)
    direction = (q_rand - q_near) / norm(q_rand - q_near);
    q_new = q_near + diag(delta_q) * direction;
end

