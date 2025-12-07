import random
import math

def generate_galaxy(filename, num_bodies=500):
    with open(filename, 'w') as f:
        f.write(f"{num_bodies}\n")
        
        # We will create two galaxies
        # Galaxy 1: Centered at (1.0, 2.0), moving right
        # Galaxy 2: Centered at (3.0, 2.0), moving left
        
        g1_center = (1.2, 2.0)
        g2_center = (2.8, 2.0)
        
        for i in range(num_bodies):
            # Split bodies between the two galaxies
            if i < num_bodies // 2:
                cx, cy = g1_center
                vx_bulk, vy_bulk = 0.5, 0.2  # Moving Right and Up
                is_g1 = True
            else:
                cx, cy = g2_center
                vx_bulk, vy_bulk = -0.5, -0.2 # Moving Left and Down
                is_g1 = False
            
            # Random position within a radius
            angle = random.uniform(0, 2 * math.pi)
            dist = random.uniform(0.05, 0.6) # Radius 0.6
            
            x = cx + dist * math.cos(angle)
            y = cy + dist * math.sin(angle)
            
            # Orbital Velocity Calculation (v = sqrt(GM/r))
            # We treat the center as having a "virtual" heavy mass to encourage rotation
            # G=1.0 in your sim. Virtual Mass ~ 5.0 (tweak for stability)
            virtual_mass = 2.0
            velocity_magnitude = math.sqrt(virtual_mass / dist)
            
            # Determine orbital direction (tangent to radius)
            # (-sin, cos) is counter-clockwise
            vx_orb = -math.sin(angle) * velocity_magnitude
            vy_orb = math.cos(angle) * velocity_magnitude
            
            # Flip rotation for second galaxy for chaos
            if not is_g1:
                vx_orb *= -1
                vy_orb *= -1
            
            # Combine Bulk + Orbital velocity
            vx = vx_bulk + vx_orb
            vy = vy_bulk + vy_orb
            
            # Mass (small particles)
            mass = 0.01
            
            # Heavier core particles for stability (optional)
            if i == 0 or i == num_bodies // 2:
                x, y = cx, cy
                vx, vy = vx_bulk, vy_bulk
                mass = 10.0 # Heavy center
            
            # Write to file
            # Index x y mass vx vy
            f.write(f"{i} {x:.6f} {y:.6f} {mass:.6f} {vx:.6f} {vy:.6f}\n")

if __name__ == "__main__":
    generate_galaxy("in/galaxy.txt", 1000)
    print("Generated in/galaxy.txt with 1000 bodies.")