include("./tree.jl")

type BBFMM
    tree::Tree
    
    chargeTree
    kernel
    rank
    
    R # translation from cluster to parent
    
    
    nCNode
    cNode # std cnode
    TNode # chey poly at cnodes
    
    BBFMM() = new(Tree(), [], Nil, 0, [], [],[],[])
end


function initialize!(b::BBFMM,  nCNode, source, target, charge, nSource, nTarget, rank, maxLevel, kernel::Function)
    b.tree = Tree()
    t = b.tree
    populate!(t, source, target, nSource, nTarget, rank, maxLevel)
    
    b.nCNode = nCNode
    b.chargeTree = charge
    b.kernel = kernel    
    b.cNode = getStandardCNodes(nCNode)    
    b.TNode = getStandardCPoly(nCNode, nCNode, b.cNode)
    b.R = getTransfer(nCNode, b.cNode, b.TNode)
    b.rank = nCNode^2


end

function FMM!(b::BBFMM)
    upPass!(b, 1)
    potential = zeros(b.tree.nTarget)
    downPass!(b, 1, potential)
    return potential
end

function upPass!(b::BBFMM, rootId)
    t = b.tree
    d = t.dict
    node = d[rootId]
    
    """
    allocate space for each node's properties.
    
    scaledCnode # cnodes
    cnodeCharge # charge on cnodes
    cnodePotential # potential on cnodes
    """
    node.cnodePotential = zeros(b.rank)
    node.cnodeCharge = zeros(b.rank)
    
    node.scaledCnode = []
    getScaledCNodes(b.nCNode, b.cNode, node.center, node.radius,node.scaledCnode)
    
    if (node.isLeaf)
        getCharge!(b, rootId)
        node.R = getTransferP2C(b.nCNode, node.source, node.center, node.radius, b.cNode, b.TNode)
        node.L = getTransferP2C(b.nCNode, node.target, node.center, node.radius, b.cNode, b.TNode)
        

        """
        interpolation 
        """
        node.cnodeCharge += node.R' * node.charge

    else
        for l = 1:4
            upPass!(b, node.child[l])
            if (!d[node.child[l]].isEmpty)
                node.cnodeCharge += (b.R[l])' * d[node.child[l]].cnodeCharge
            end
        end
    end
end


function downPass!(b::BBFMM, rootId, potential)
    t = b.tree
    d = t.dict
    node = d[rootId]
    
    if (node.parent != -1)
        """
        Vlist
        """
        for l = 1:node.nVList
            if (!d[node.vList[l]].isEmpty)
                node.cnodePotential += getAccumulationCheb(b, b.nCNode, d[node.vList[l]].scaledCnode, b.nCNode, node.scaledCnode ,d[node.vList[l]].cnodeCharge)
            end
        end
        """
        XList
        """

        for l = 1:node.nXList
            if (!d[node.xList[l]].isEmpty)
                node.cnodePotential += getAccumulationCheb(b, b.nCNode, d[node.xList[l]].scaledCnode, b.nCNode, node.scaledCnode,d[node.xList[l]].cnodeCharge)
            end
        end
        """
        L2L
        """
        parentNode = d[node.parent]
        node.cnodePotential += b.R[node.nodeIndex]*parentNode.cnodePotential
       
    end
    
    if (node.isLeaf && node.nTarget!=0)
        node.potential = zeros(node.nTarget)
        
        """
        UList
        """
        for l = 1:node.nUList
            if(!d[node.uList[l]].isEmpty)
                getCharge!(b, node.uList[l])
                node.potential += getAccumulation(b, d[node.uList[l]].source, node.target, d[node.uList[l]].charge)
            end
        end
        
        """
        WList, should not be direct
        """
        
        for l = 1:node.nWList
            if(!d[node.wList[l]].isEmpty)
                getCharge!(b, node.wList[l])
                node.potential += getAccumulation(b, d[node.wList[l]].source, node.target, d[node.wList[l]].charge)
            end
        end
        
        """
        from potential on cheb nodes
        """
        node.potential += node.L * node.cnodePotential
        
        """
        finalize
        """
        for l = 1:node.nTarget
            potential[node.targetIndex[l]] += node.potential[l]
        end
    end
    
    if !node.isLeaf
        for l = 1:4
            downPass!(b, node.child[l], potential)
        end
    end
end


function downPass2!(b::BBFMM, rootId, potential)
    t = b.tree
    d = t.dict
    node = d[rootId]
    if (!node.isEmpty)
        if (node.isLeaf)
            #ULIST
            node.potential = zeros(node.nTarget)
            for l = 1:node.nUList
                if(!d[node.uList[l]].isEmpty)
                    getCharge!(b, node.uList[l])
                    node.potential += getAccumulation(b, d[node.uList[l]].source, node.target, d[node.uList[l]].charge)
                end
            end
            # w empty
            
            node.potential += node.R*node.cnodePotential
            for l = 1:node.nTarget
                potential[node.targetIndex[l]] += node.potential[l]
            end
        else        

            for k = 1:4
                childId = node.child[k]
                childnode = d[childId]
                if (!childnode.isEmpty)
                    for l = 1:childnode.nVList
                        childnode.cnodePotential += getAccumulationCheb(b, b.nCNode, d[childnode.vList[l]].scaledCnode, b.nCNode, childnode.scaledCnode ,d[childnode.vList[l]].cnodeCharge)
 
                    end
                end
            end
  
            for k = 1:4
                childId = node.child[k]
                childnode = d[childId]
                if (!childnode.isEmpty)
                    childnode.cnodePotential += b.R[k] * node.cnodePotential
                end
            end
            
            for k = 1:4
                downPass2!(b, node.child[k], potential)
            end
        end
    end
end

    
function getCharge!(b::BBFMM, id)
    # lazy evaluation
    node = b.tree.dict[id]
    if (node.chargeComputed == true)
        # pass
    else
        node.chargeComputed = true
        node.charge = zeros(node.nSource)
        for l = 1:node.nSource
            node.charge[l] = b.chargeTree[node.sourceIndex[l]]
        end
    end
end


function getScaledCNodes(nCNodes, cNode, center, radius, cheNode)
    for k = 1:nCNodes
        push!(cheNode, center + cNode[k]*radius)
    end
end


function getStandardCNodes(nCNodes)
    cNode = zeros(nCNodes)
    for k = 1:nCNodes
        cNode[k] = -cos((k - 0.5)*π/nCNodes)
    end
    return cNode
end

function getStandardCPoly(nPoly, N, x)
    Tnode = zeros(N, nPoly)
    Tnode[:, 1] = ones(N)
    
    if (nPoly > 1)
        Tnode[:, 2] = x
        for k = 3:nPoly 
            Tnode[:, k] = 2.0 * x.*Tnode[:, k-1] - Tnode[:, k-2]
        end
    end
    
    return Tnode
end
   

function getTransfer(nCNode, cNode, TNode)
    S = getTransferP2C_CNode(nCNode, cNode, TNode)
    Transfer0 = S[1:nCNode, 1:nCNode]
    Transfer1 = S[(nCNode+1):2*nCNode, 1:nCNode]
    rank = nCNode * nCNode
    R1 = zeros(rank, rank)
    R2 = zeros(rank, rank)
    R3 = zeros(rank, rank)
    R4 = zeros(rank, rank)
    for i = 1:nCNode
        for j = 1:nCNode
            for k = 1:nCNode
                for l=1:nCNode
                    R1[(i-1)*nCNode+j, (k-1)*nCNode+l] = Transfer0[i,k]*Transfer0[j,l]
                    R2[(i-1)*nCNode+j, (k-1)*nCNode+l] = Transfer0[i,k]*Transfer1[j,l]
                    R3[(i-1)*nCNode+j, (k-1)*nCNode+l] = Transfer1[i,k]*Transfer0[j,l]
                    R4[(i-1)*nCNode+j, (k-1)*nCNode+l] = Transfer1[i,k]*Transfer1[j,l]
                end
            end
        end
    end
    return [R1, R2, R3, R4]
end

function getTransferP2C_CNode(nCNode, cNode, TNode)
    childcNode = zeros(2 * nCNode)
    childcNode[1:nCNode] = 0.5 * (cNode - ones(nCNode))
    childcNode[(nCNode + 1): 2*nCNode] = 0.5 * (cNode + ones(nCNode))
    Transfer = getStandardCPoly(nCNode, 2 * nCNode, childcNode)
    return (2.0 * Transfer * TNode' - ones(2*nCNode, nCNode))/nCNode
end

function getTransferP2C(nCNode, source, center, radius, cNode, TNode)
    N = length(source)
    loc0 = zeros(N)
    loc1 = zeros(N)
    for i = 1:N
        loc0[i] = (source[i][1] - center[1])/radius[1]
        loc1[i] = (source[i][2] - center[2])/radius[2]
    end
    
    T0 = getStandardCPoly(nCNode, N, loc0)
    T1 = getStandardCPoly(nCNode, N, loc1)
 
    T0 = (2.0 * T0 * TNode' - ones(N, nCNode))/nCNode
    T1 = (2.0 * T1 * TNode' - ones(N, nCNode))/nCNode
        
    rank = nCNode^2
    R = zeros(N, rank)
    for k = 1:N
        for i = 1:nCNode
            for j = 1:nCNode
                R[k, (j-1)*nCNode + i] = T0[k,i]*T1[k,j]
            end
        end
    end
    return R
end


function kernelEval(source, target, f)
    sourceSize = length(source)
    targetSize = length(target)
    K = zeros(targetSize, sourceSize)
    for s = 1:sourceSize
        for t = 1: targetSize
            K[t, s] = f(source[s], target[t])
        end
    end
    return K
end


function getAccumulation(b::BBFMM, source, target, charge)
    K = kernelEval(source, target, b.kernel)
    # K(target, source) * charge
    # density of length(target)
    # charge of lenth(source)
    return K * charge
end

function getAccumulationCheb(b::BBFMM, M, s, N, t, charge)
    sourceSize = length(s)^2
    targetSize = length(t)^2
    source = []
    target = []
    for j = 1:M
        for i = 1:M
            push!(source, [s[i][1], s[j][2]])
        end
    end
    
    for j = 1:N
        for i=1:N
            push!(target, [t[i][1], t[j][2]])
        end
    end
    
    K = kernelEval(source, target, b.kernel)
    return K*charge

end


