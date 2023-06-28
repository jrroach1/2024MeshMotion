"""
    cylindermesh(; ref, porder, bunch=2.5, format="gmsh")

Create high-order cylinder mesh for the Moving Mesh test case in the CFD Convergence Workshop.

`ref`:     uniform refinement index (0, 1, 2, 3, ...) 

`porder`:  curved mesh order (1 = linear, 2, 3, ...)

`bunch`:   BL bunching factor (2.5 default, higher -> more bunching on cylinder wall)

# Examples
```
# Generate file cyl_ref1_p3_b250.msh (from within Julia):
include("cylindermesh.jl")
cylindermesh(ref=1, porder=3);

# Same, from command-line:
julia -L cylindermesh.jl -e 'cylindermesh(ref=1, porder=3)'

# From command-line, use gmsh to convert to latest format:
gmsh -convert cyl_ref1_p3_b250.msh

# From command-line, use gmsh to convert to CGNS format:
gmsh -format cgns -save cyl_ref1_p3_b250.msh
```
"""
function cylindermesh(; ref, porder, bunch=2.5, format="gmsh")
    
    # GMSH format
    gmsh_lin_types = [1,8,26,27,28]
    gmsh_quad_types = [3,10,36,37,38]
    gmsh_lin_map = [ [1,2],
                     [1,3,2],
                     [1,4,2,3],
                     [1,5,2,3,4],
                     [1,6,2,3,4,5] ]

    gmsh_quad_map = [
    [1, 2, 4, 3],
    [1, 5, 2, 8, 9, 6, 4, 7, 3],
    [1, 5, 6, 2, 12, 13, 14, 7, 11, 16, 15, 8, 4, 10, 9, 3],
    [1, 5, 6, 7, 2, 16, 17, 21, 18, 8, 15, 24, 25, 22, 9, 14, 20, 23, 19, 10, 4, 13, 12, 11, 3],
    [1, 5, 6, 7, 8, 2, 20, 21, 25, 26, 22, 9, 19, 32, 33, 34, 27, 10, 18, 31, 36, 35, 28, 11, 17, 24, 30, 29, 23, 12, 4, 16, 15, 14, 13, 3] ]

    # Utilities

    function gindex(NBG,g,n)
        if g<0
            g = -g
            n = size(NBG,2) + 1 .- n
        end
        return I = NBG[g,n]
    end

    function mesh_edge(V,N,ne)
        s = range(0, 1, length=ne+1)
        nV = size(V,1)
        v = zeros(ne+1,2)
        for i = 1:2
            v[:,i] = V[N[1],i] * (1 .- s) + V[N[2],i]*s
        end
        I = [N[1]; nV+1:nV+ne-1; N[2]]'
        V = [V; v[2:ne,:]]
        return V,I
    end
    
    function mesh_elem(V,NBG,G)
        nn = size(NBG,2)
        ne = nn-1
        I = zeros(Int64, 1, nn*nn)
        I[1:nn] = gindex(NBG,G[1],1:nn)
        for j = 2:nn-1
            V,J = mesh_edge(V, [gindex(NBG,G[3],j); gindex(NBG,G[4],j)], ne)
            I[(j-1)*nn .+ (1:nn)] = J
        end
        I[end-nn+1:end] = gindex(NBG,G[2], 1:nn)
        return V,I
    end

    # Sizes
    ne0 = 2^(ref+1)  # number of intervals per edge
    ne = ne0*porder   # including high-order nodes
    
    # Geometric parameters
    l = 1 # half-edge length
    f = 0.5 # scaling factor
    R0 = l*sqrt(2) # diagonal
    R = 0.5  # cylinder radius
    
    # define blocks
    BV = [-l -l; l -l; -l l; l l]
    BV = [BV; BV*f] # block vertices
    BG = [1 2; 1 3; 2 4; 3 4; 5 6; 5 7; 6 8; 7 8; 1 5; 2 6; 3 7; 4 8] # block edges
    BE = [1 5 9 10; 9 11 2 6; -10 -12 7 3; 8 4 -11 -12; 5 8 6 7] # block elements
    BB = [1 2 3 4] # block boundaries

    # global vertex list, start with block vertices
    V = BV
    nV = size(V,1)

    # mesh block edges
    nbg = size(BG,1)
    NBG = zeros(Int64, nbg, ne+1)
    for g = 1:nbg
        V,I = mesh_edge(V,BG[g,:],ne)
        NBG[g,:] = I
    end

    # mesh block elements
    nbe = size(BE,1)
    NBE = zeros(Int64, nbe, (ne+1)^2)
    for e = 1:nbe
        V,I = mesh_elem(V,NBG,BE[e,:])
        NBE[e,:] = I
    end

    # transform to cylinder shape
    for i=1:size(V,1)
        r = sqrt(sum(V[i,:].^2))
        t0 = atan(V[i,2],V[i,1])
        t = t0 - pi/2*floor(t0*2/pi)
        t1 = t
        t = min(t, pi/2-t)
        d = l/cos(t)
        f = r/d
        dd = R0 - d
        r = r + dd*f^2*sin(f*pi/2)
        b = r / R0
        a = bunch
        r = R0 * (1-exp(-a*b)) / (1-exp(-a))
        t = t0 - 0.06*f*sin(t1*4)
        V[i,:] = r * [cos(t), sin(t)] * R/R0
    end
    
    nelem = nbe * ne0^2
    Qp = porder+1
    nV = size(V,1)

    fname = "cyl_ref$(ref)_p$(porder)_b$(round(Int,100bunch))"

    if format == "gmsh"
        fid = open(fname * ".msh", "w")

        println(fid, """
            \$MeshFormat
            2.2 0 8
            \$EndMeshFormat
            \$PhysicalNames
            2
            1 1 "Cylinder Boundary"
            2 2 "Interior"
            \$EndPhysicalNames
            \$Nodes""")

        println(fid, "$nV")
        for i = 1:nV
            println(fid,"$i $(V[i,1]) $(V[i,2]) 0")
        end

        println(fid, """
             \$EndNodes
             \$Elements""")

        nbnd = ne0*length(BB)
        println(fid, "$(nbnd+nelem)")
        el = 0

        el_type = gmsh_lin_types[porder]
        for b = 1:length(BB)
            for g = 1:porder:ne
                el += 1
                print(fid, "$el $el_type 2 1 1")
                n = NBG[BB[b],g:g+porder]
                n .= n[gmsh_lin_map[porder]]
                for nn in n
                    print(fid, " $nn")
                end
                println(fid)
            end
        end
        
        el_type = gmsh_quad_types[porder]
        for e=1:nbe
            for j=1:ne0
                for i=1:ne0
                    el += 1
                    print(fid, "$el $el_type 2 2 2")
                    n0 = (j-1)*(ne+1)*porder + (i-1)*porder
                    n = zeros(Int64, Qp*Qp)
                    k = 0
                    for jq=1:Qp
                        for iq=1:Qp
                            k = k+1
                            n[k] = n0 + (jq-1)*(ne+1) + iq
                        end
                    end
                    n[gmsh_quad_map[porder]] .= n
                    for k=1:length(n)
                        print(fid, " $(NBE[e,n[k]])")
                    end
                    println(fid)
                end
            end
        end
        println(fid, "\$EndElements")
        
        close(fid)
    elseif format == "gri"
        fid = open(fname * ".gri", "w")
        println(fid, "$(size(V,1)) $nelem 2")
        for i = 1:size(V,1)
            println(fid, "$(V[i,1]) $(V[i,2])")
        end
        println(fid, "1") # 1 boundary group
        println(fid, "$(ne0*length(BB)) 2 Cylinder") # bfaces per boundary group
        for b = 1:length(BB)
            for g = 1:porder:ne
                println(fid, "$(NBG[BB[b],g]) $(NBG[BB[b],g+porder])")
            end
        end
        
        println(fid, "$nelem $porder QuadLagrange")
        for e=1:nbe
            for j=1:ne0
                for i=1:ne0
                    n0 = (j-1)*(ne+1)*porder + (i-1)*porder
                    n = zeros(Int64, Qp*Qp)
                    k = 0
                    for jq=1:Qp
                        for iq=1:Qp
                            k = k+1
                            n[k] = n0 + (jq-1)*(ne+1) + iq
                        end
                    end
                    for k=1:length(n)
                        print(fid, "$(NBE[e,n[k]]) ")
                    end
                    println(fid)
                end
            end
        end
        
        close(fid)
    else
        error("Unknown format")
    end
    
    return V,NBE,NBG
end
