"""
rotinas_talude.py
=================
Rotinas complementares a `rotinas_MEF.py` para a analise de tensoes e
deformacoes de um talude sob peso proprio (elementos LST, deformacao plana).

Conteudo
--------
Geometria / malha
    malha_talude_lst        malha LST estruturada (mapeamento vertical) para
                            um talude definido pela linha da superficie
Carregamento e solucao
    forcas_peso_proprio     vetor de forcas nodais equivalentes ao peso proprio
    montar_K_esparsa        montagem esparsa da matriz de rigidez (usa k_lst)
    SolverApoios            fatoracao LU esparsa com graus de liberdade fixos
    montar_B_gauss          matrizes B em todos os pontos de Gauss (usa b_lst)
Plasticidade
    retorno_mc              Mohr-Coulomb: retorno das tensoes (plano principal,
                            arestas e apice), fluxo nao associado, vetorizado
    mobilizacao_mc          grau de mobilizacao da resistencia
Pos-processamento
    invariantes             tensoes principais, tau_max, p e mobilizacao
    plot_campo_gauss        campo escalar por pontos de Gauss (sub-areas)
    plot_campo_nodal        campo nodal (contornos suaves)
    plot_direcoes_principais  setas das direcoes principais de tensao

Convencoes (iguais as de `rotinas_MEF.py`)
    * No LST: 1-3 vertices (anti-horario), 4 entre 1-2, 5 entre 2-3, 6 entre 3-1
    * Graus de liberdade do elemento: [u1 v1 u2 v2 ... u6 v6]
    * Vetor de tensoes (deformacao plana): [sx, sy, sz, txy]
    * Tensoes: TRACAO POSITIVA, COMPRESSAO NEGATIVA (convencao da Mecanica)
    * Unidades sugeridas: kN, m, kPa, kN/m3
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Polygon
import matplotlib.tri as mtri
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import splu

from rotinas_MEF import b_lst, k_lst, gauss_tri, gauss_akoor


# =============================================================================
#  GEOMETRIA / MALHA
# =============================================================================
def malha_talude_lst(superficie, y_base, dx_alvo, ny=None):
    """
    Gera uma malha estruturada de elementos LST (6 nos) para um talude.

    A malha e obtida por mapeamento vertical: a linha da superficie e dividida
    em colunas (com nos nos vertices informados) e cada coluna e dividida em
    `ny` camadas entre y_base e a superficie. Cada quadrilatero e dividido em
    dois triangulos (diagonais alternadas) e os nos do meio dos lados sao
    criados em seguida.

    Parametros
    ----------
    superficie : array_like (n, 2)
        Pontos (x, y) da linha da superficie, com x estritamente crescente
        (nao sao admitidos trechos verticais).
    y_base : float
        Cota da base (fronteira inferior do dominio).
    dx_alvo : float
        Tamanho alvo (em x) das colunas de elementos [m].
    ny : int, opcional
        Numero de camadas. Se omitido: ceil(altura_maxima / dx_alvo).

    Retorna
    -------
    ElemNode : ndarray (nelem, 6) de int, numeracao 1-based (como em
               `lst_topo_taper_beam`)
    SysCoor  : ndarray (nnode, 2)
    NodeSets : dict (1-based) com os conjuntos "left", "right", "base",
               "surface" (todos os nos, inclusive os do meio dos lados)
    info     : dict com "nx", "ny", "x_col" (abscissas das colunas) e
               "y_surf_col" (cota da superficie em cada coluna)
    """
    sup = np.asarray(superficie, dtype=float)
    if sup.ndim != 2 or sup.shape[1] != 2 or sup.shape[0] < 2:
        raise ValueError("`superficie` deve ter ao menos 2 pontos (x, y).")
    if np.any(np.diff(sup[:, 0]) <= 0):
        raise ValueError("A abscissa x da superficie deve ser estritamente crescente.")
    if np.any(sup[:, 1] <= y_base):
        raise ValueError("Todos os pontos da superficie devem estar acima de y_base.")

    # ----- colunas: subdivide cada trecho da superficie -------------------
    x_col = [sup[0, 0]]
    for k in range(sup.shape[0] - 1):
        seg = sup[k + 1, 0] - sup[k, 0]
        n = max(1, int(np.ceil(seg / dx_alvo - 1e-9)))
        x_col.extend(list(sup[k, 0] + seg * np.arange(1, n + 1) / n))
    x_col = np.asarray(x_col)
    y_surf = np.interp(x_col, sup[:, 0], sup[:, 1])
    nx = x_col.size - 1

    if ny is None:
        ny = max(2, int(np.ceil((y_surf.max() - y_base) / dx_alvo)))

    # ----- nos de vertice -------------------------------------------------
    nodes = []
    for j in range(ny + 1):                      # camadas (de baixo p/ cima)
        for i in range(nx + 1):
            y = y_base + (y_surf[i] - y_base) * j / ny
            nodes.append((x_col[i], y))
    nid = lambda i, j: j * (nx + 1) + i          # indice (0-based) do no (i, j)

    # ----- triangulos (anti-horarios) ------------------------------------
    tris = []
    for j in range(ny):
        for i in range(nx):
            n00, n10 = nid(i, j), nid(i + 1, j)
            n11, n01 = nid(i + 1, j + 1), nid(i, j + 1)
            if (i + j) % 2 == 0:
                tris.append((n00, n10, n11))
                tris.append((n00, n11, n01))
            else:
                tris.append((n00, n10, n01))
                tris.append((n10, n11, n01))

    # ----- nos do meio dos lados -----------------------------------------
    nodes = list(nodes)
    edge_mid = {}

    def meio(a, b):
        key = (a, b) if a < b else (b, a)
        if key not in edge_mid:
            xa, ya = nodes[a]
            xb, yb = nodes[b]
            nodes.append((0.5 * (xa + xb), 0.5 * (ya + yb)))
            edge_mid[key] = len(nodes) - 1
        return edge_mid[key]

    elem = []
    for (a, b, c) in tris:
        elem.append((a, b, c, meio(a, b), meio(b, c), meio(c, a)))

    ElemNode0 = np.asarray(elem, dtype=int)
    SysCoor = np.asarray(nodes, dtype=float)

    # sanidade: area positiva (anti-horario) em todos os elementos
    X = SysCoor[ElemNode0[:, :3], 0]
    Y = SysCoor[ElemNode0[:, :3], 1]
    area = 0.5 * ((X[:, 1] - X[:, 0]) * (Y[:, 2] - Y[:, 0])
                  - (X[:, 2] - X[:, 0]) * (Y[:, 1] - Y[:, 0]))
    if np.any(area <= 0):
        raise RuntimeError("Elementos com area nao positiva: revise a geometria.")

    # ----- conjuntos de nos (1-based) -------------------------------------
    tol = 1e-9 * max(1.0, np.ptp(SysCoor[:, 0]))
    ysurf_nodes = np.interp(SysCoor[:, 0], sup[:, 0], sup[:, 1])
    NodeSets = {
        "left": np.where(np.abs(SysCoor[:, 0] - x_col[0]) < tol)[0] + 1,
        "right": np.where(np.abs(SysCoor[:, 0] - x_col[-1]) < tol)[0] + 1,
        "base": np.where(np.abs(SysCoor[:, 1] - y_base) < tol)[0] + 1,
        "surface": np.where(np.abs(SysCoor[:, 1] - ysurf_nodes) < tol)[0] + 1,
    }
    info = {"nx": nx, "ny": ny, "x_col": x_col, "y_surf_col": y_surf}
    return ElemNode0 + 1, SysCoor, NodeSets, info


# =============================================================================
#  CARREGAMENTO, MONTAGEM E SOLUCAO
# =============================================================================
def forcas_peso_proprio(ElemX, ElemY, ElemDof, ndof, gamma, t, gaussorder):
    """
    Forcas nodais equivalentes ao peso proprio (forca de corpo b = (0, -gamma)):

        f_e = integral( N^T b ) dV = sum_g W_g * N(xi_g)^T * b * A * t

    Retorna o vetor global f (ndof,). Para um elemento LST os nos de canto
    recebem parcela nula e cada no do meio 1/3 do peso do elemento (integracao
    exata com 3 ou mais pontos de Gauss).
    """
    XI, W = gauss_tri(gaussorder)
    W = np.asarray(W).reshape(-1)
    x1, x2, x3 = XI[:, 0], XI[:, 1], XI[:, 2]
    # Funcoes de forma do LST nos pontos de Gauss: (ng, 6)
    N = np.column_stack((x1 * (2 * x1 - 1), x2 * (2 * x2 - 1), x3 * (2 * x3 - 1),
                         4 * x1 * x2, 4 * x2 * x3, 4 * x3 * x1))
    area = 0.5 * ((ElemX[:, 1] - ElemX[:, 0]) * (ElemY[:, 2] - ElemY[:, 0])
                  - (ElemX[:, 2] - ElemX[:, 0]) * (ElemY[:, 1] - ElemY[:, 0]))
    Nw = W @ N                                   # (6,) = integral das N / A
    fy = -gamma * t * area[:, None] * Nw[None, :]            # (nelem, 6)
    f = np.zeros(ndof)
    np.add.at(f, ElemDof[:, 1::2], fy)           # graus de liberdade v (y)
    return f


def montar_K_esparsa(ElemX, ElemY, ElemDof, D, t, gaussorder, s_type):
    """Monta a matriz de rigidez global (esparsa, CSC) com `k_lst`."""
    nelem = ElemX.shape[0]
    ndof = ElemDof.max() + 1
    nloc = ElemDof.shape[1]
    rows = np.empty(nelem * nloc * nloc, dtype=int)
    cols = np.empty_like(rows)
    vals = np.empty(nelem * nloc * nloc)
    p = 0
    for e in range(nelem):
        Ke = k_lst(np.column_stack((ElemX[e], ElemY[e])), D, t, gaussorder, s_type)
        d = ElemDof[e]
        rows[p:p + nloc * nloc] = np.repeat(d, nloc)
        cols[p:p + nloc * nloc] = np.tile(d, nloc)
        vals[p:p + nloc * nloc] = Ke.ravel()
        p += nloc * nloc
    return coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsc()


class SolverApoios:
    """
    Resolve K u = f com graus de liberdade restritos a ZERO (apoios).

    A fatoracao LU esparsa dos graus livres e feita uma unica vez e pode ser
    reutilizada em todas as iteracoes (metodo de Newton-Raphson modificado).
    As numeracoes dos graus de liberdade sao 0-based.
    """

    def __init__(self, K, dof_fixos):
        ndof = K.shape[0]
        self.ndof = ndof
        self.fixos = np.unique(np.asarray(dof_fixos, dtype=int))
        mask = np.ones(ndof, dtype=bool)
        mask[self.fixos] = False
        self.livres = np.where(mask)[0]
        self.lu = splu(csc_matrix(K[self.livres][:, self.livres]))

    def resolver(self, f):
        u = np.zeros(self.ndof)
        u[self.livres] = self.lu.solve(np.asarray(f)[self.livres])
        return u


def montar_B_gauss(ElemX, ElemY, gaussorder, s_type):
    """
    Matrizes B (b_lst) em todos os pontos de Gauss.

    Retorna
    -------
    Bg : (nelem, ng, nsigma, 12)
    wg : (nelem, ng)  pesos de integracao  W_g * A_e  (multiplicar por t)
    Xg, Yg : (nelem, ng) coordenadas dos pontos de Gauss
    """
    XI, W = gauss_tri(gaussorder)
    W = np.asarray(W).reshape(-1)
    nelem = ElemX.shape[0]
    nsigma = 3 if s_type == 1 else 4
    Bg = np.zeros((nelem, gaussorder, nsigma, 12))
    for e in range(nelem):
        coor = np.column_stack((ElemX[e], ElemY[e]))
        for g in range(gaussorder):
            Bg[e, g] = b_lst(coor, XI[g], s_type)
    area = 0.5 * ((ElemX[:, 1] - ElemX[:, 0]) * (ElemY[:, 2] - ElemY[:, 0])
                  - (ElemX[:, 2] - ElemX[:, 0]) * (ElemY[:, 1] - ElemY[:, 0]))
    wg = area[:, None] * W[None, :]
    Xg = ElemX[:, :3] @ XI.T
    Yg = ElemY[:, :3] @ XI.T
    return Bg, wg, Xg, Yg


# =============================================================================
#  PLASTICIDADE (Mohr-Coulomb, retorno no espaco das tensoes principais)
# =============================================================================
def _principais(sig):
    """
    Tensoes principais de [sx, sy, sz, txy] (deformacao plana: sz e principal).

    Retorna P (..., 3) = [p1, p2, sz] (p1 >= p2 no plano), o centro e o raio do
    circulo de Mohr no plano e cos(2 theta), sin(2 theta) da direcao de p1.
    """
    sx, sy, sz, txy = sig[..., 0], sig[..., 1], sig[..., 2], sig[..., 3]
    centro = 0.5 * (sx + sy)
    raio = np.hypot(0.5 * (sx - sy), txy)
    seguro = np.where(raio > 1e-14, raio, 1.0)
    c2 = np.where(raio > 1e-14, 0.5 * (sx - sy) / seguro, 1.0)
    s2 = np.where(raio > 1e-14, txy / seguro, 0.0)
    P = np.stack((centro + raio, centro - raio, sz), axis=-1)
    return P, centro, raio, c2, s2


def retorno_mc(sig_tr, E, nu, c, phi_graus, psi_graus=0.0):
    """
    Retorno das tensoes ao criterio de Mohr-Coulomb (vetorizado), no espaco das
    tensoes principais, com lei de fluxo nao associada (dilatancia psi):

        f = (1 + sin phi) s1 - (1 - sin phi) s3 - 2 c cos phi <= 0
        g = (1 + sin psi) s1 - (1 - sin psi) s3         (s1 >= s2 >= s3, tracao +)

    Trata plano principal, arestas (duas superficies ativas) e apice do cone
    (de Souza Neto, Peric e Owen, 2008, cap. 8; Clausen et al., 2007).

    Parametros
    ----------
    sig_tr : (..., 4) tensoes tentativa (elasticas) [sx, sy, sz, txy]

    Retorna
    -------
    sig : (..., 4) tensoes atualizadas
    plastico : (...) bool, pontos que escoaram nesta atualizacao
    """
    sphi = np.sin(np.radians(phi_graus))
    cphi = np.cos(np.radians(phi_graus))
    spsi = np.sin(np.radians(min(psi_graus, phi_graus)))
    G = E / (2.0 * (1.0 + nu))
    Kb = E / (3.0 * (1.0 - 2.0 * nu))
    lamL = Kb - 2.0 * G / 3.0
    coh = 2.0 * c * cphi

    P, centro, raio, c2, s2 = _principais(sig_tr)
    ordem = np.argsort(-P, axis=-1)                       # decrescente
    S = np.take_along_axis(P, ordem, axis=-1)
    s1, s2_, s3 = S[..., 0], S[..., 1], S[..., 2]

    f_tr = (1 + sphi) * s1 - (1 - sphi) * s3 - coh
    plastico = f_tr > 0.0
    if not np.any(plastico):
        return sig_tr.copy(), plastico

    # --- plano principal ---------------------------------------------------------------
    Maa = 4 * G * (1 + sphi * spsi) + 4 * lamL * sphi * spsi
    dl = f_tr / Maa
    A_ = 2 * G * (1 + spsi) + 2 * lamL * spsi              # reducao de s1
    C_ = 2 * G * (1 - spsi) - 2 * lamL * spsi              # aumento de s3
    n1 = s1 - dl * A_
    n2 = s2_ - 2 * lamL * spsi * dl
    n3 = s3 + dl * C_
    cruza_dir = n1 < n2                                    # s1 cruzou s2 -> aresta s1 = s2
    cruza_esq = n3 > n2                                    # s3 cruzou s2 -> aresta s2 = s3

    # --- aresta direita (s1 = s2) ----------------------------------------------------------
    Mab_d = 2 * G * (1 - sphi) * (1 - spsi) + 4 * lamL * sphi * spsi
    fa = f_tr
    fb = (1 + sphi) * s2_ - (1 - sphi) * s3 - coh
    det_d = Maa ** 2 - Mab_d ** 2
    la_d = (Maa * fa - Mab_d * fb) / det_d
    lb_d = (Maa * fb - Mab_d * fa) / det_d
    vol_d = 2 * lamL * spsi * (la_d + lb_d)
    d1 = s1 - 2 * G * (1 + spsi) * la_d - vol_d
    d2 = s2_ - 2 * G * (1 + spsi) * lb_d - vol_d
    d3 = s3 + 2 * G * (1 - spsi) * (la_d + lb_d) - vol_d
    ok_d = (la_d >= 0) & (lb_d >= 0) & (d1 >= d3)          # multiplicadores >= 0 e ordem valida

    # --- aresta esquerda (s2 = s3) ----------------------------------------------------------
    Mac_e = 2 * G * (1 + sphi) * (1 + spsi) + 4 * lamL * sphi * spsi
    fc = (1 + sphi) * s1 - (1 - sphi) * s2_ - coh
    det_e = Maa ** 2 - Mac_e ** 2
    la_e = (Maa * fa - Mac_e * fc) / det_e
    lc_e = (Maa * fc - Mac_e * fa) / det_e
    vol_e = 2 * lamL * spsi * (la_e + lc_e)
    e1 = s1 - 2 * G * (1 + spsi) * (la_e + lc_e) - vol_e
    e2 = s2_ + 2 * G * (1 - spsi) * lc_e - vol_e
    e3 = s3 + 2 * G * (1 - spsi) * la_e - vol_e
    ok_e = (la_e >= 0) & (lc_e >= 0) & (e1 >= e3)

    # --- apice (tensao isotropica c / tan phi) -------------------------------------------------
    if sphi > 1e-12:
        p_ap = c * cphi / sphi
    else:                                                   # phi = 0: sem apice
        p_ap = np.inf

    # --- escolha da regiao de retorno ---------------------------------------------------------
    usa_dir = plastico & cruza_dir
    usa_esq = plastico & ~cruza_dir & cruza_esq
    usa_pl = plastico & ~cruza_dir & ~cruza_esq
    usa_ap = (usa_dir & ~ok_d) | (usa_esq & ~ok_e)
    usa_dir &= ~usa_ap
    usa_esq &= ~usa_ap

    R1 = np.where(usa_pl, n1, np.where(usa_dir, d1, np.where(usa_esq, e1, s1)))
    R2 = np.where(usa_pl, n2, np.where(usa_dir, d2, np.where(usa_esq, e2, s2_)))
    R3 = np.where(usa_pl, n3, np.where(usa_dir, d3, np.where(usa_esq, e3, s3)))
    if np.isfinite(p_ap):
        R1 = np.where(usa_ap, p_ap, R1)
        R2 = np.where(usa_ap, p_ap, R2)
        R3 = np.where(usa_ap, p_ap, R3)
    else:                                                   # phi = 0 (Tresca): mantem media
        pm = (s1 + s2_ + s3) / 3.0
        R1 = np.where(usa_ap, pm, R1)
        R2 = np.where(usa_ap, pm, R2)
        R3 = np.where(usa_ap, pm, R3)

    # --- volta a ordem original e reconstroi [sx, sy, sz, txy] ---------------------------------
    Rs = np.stack((R1, R2, R3), axis=-1)
    Pn = np.empty_like(Rs)
    np.put_along_axis(Pn, ordem, Rs, axis=-1)
    centro_n = 0.5 * (Pn[..., 0] + Pn[..., 1])
    raio_n = 0.5 * (Pn[..., 0] - Pn[..., 1])
    sig = np.stack((centro_n + raio_n * c2, centro_n - raio_n * c2, Pn[..., 2],
                    raio_n * s2), axis=-1)
    return np.where(plastico[..., None], sig, sig_tr), plastico


def mobilizacao_mc(sig, c, phi_graus):
    """
    Grau de mobilizacao da resistencia de Mohr-Coulomb:
        ((s1 - s3)/2) / (c cos(phi) - ((s1 + s3)/2) sin(phi))
    com s1 >= s3 as tensoes principais extremas (incluindo sz). Vale 1 sobre a
    superficie de ruptura.
    """
    P, *_ = _principais(sig)
    s1, s3 = P.max(axis=-1), P.min(axis=-1)
    sphi, cphi = np.sin(np.radians(phi_graus)), np.cos(np.radians(phi_graus))
    resist = c * cphi - 0.5 * (s1 + s3) * sphi
    return np.where(resist > 1e-9, 0.5 * (s1 - s3) / np.maximum(resist, 1e-9), 1.0)


# =============================================================================
#  POS-PROCESSAMENTO
# =============================================================================
def invariantes(sig, c=None, phi_graus=None):
    """
    Grandezas derivadas das tensoes [sx, sy, sz, txy] (tracao positiva).

    Retorna dict com:
        s1, s3  tensoes principais no plano (s1 >= s3; compressao negativa)
        tmax    tensao cisalhante maxima no plano (s1 - s3)/2
        theta   angulo (graus) da direcao de s1 com o eixo x
        p       tensao media (sx+sy+sz)/3
        mob     (se c e phi_graus) grau de mobilizacao de Mohr-Coulomb
                (1 = ponto sobre a superficie de ruptura)
    """
    sx, sy, txy = sig[..., 0], sig[..., 1], sig[..., 3]
    centro = 0.5 * (sx + sy)
    raio = np.sqrt((0.5 * (sx - sy)) ** 2 + txy ** 2)
    out = {
        "s1": centro + raio,
        "s3": centro - raio,
        "tmax": raio,
        "theta": np.degrees(0.5 * np.arctan2(2 * txy, sx - sy)),
        "p": (sx + sy + sig[..., 2]) / 3.0,
    }
    if c is not None and phi_graus is not None:
        out["mob"] = mobilizacao_mc(sig, c, phi_graus)
    return out


def plot_campo_gauss(ElemX, ElemY, campo, gaussorder, ax=None, cmap="viridis",
                     vmin=None, vmax=None, titulo=None, rotulo=None,
                     malha=False):
    """
    Campo escalar por pontos de Gauss, desenhado em sub-areas (como
    `stress_plot_lst`, porem aceitando valores negativos).

    campo : (nelem, ng)
    """
    ax = plt.gca() if ax is None else ax
    Xi1, Xi2, Xi3 = gauss_akoor(gaussorder)
    nelem = ElemX.shape[0]
    polys, cores = [], []
    for e in range(nelem):
        Tx, Ty = ElemX[e, :3], ElemY[e, :3]
        for g in range(gaussorder):
            n = 3 if Xi1[g, 3] == -1 else 4
            ac = np.column_stack((Xi1[g, :n], Xi2[g, :n], Xi3[g, :n]))
            polys.append(np.column_stack((ac @ Tx, ac @ Ty)))
            cores.append(campo[e, g])
    pc = PatchCollection([Polygon(p, closed=True) for p in polys],
                         edgecolor="none", cmap=cmap)
    pc.set_array(np.asarray(cores))
    pc.set_clim(np.min(campo) if vmin is None else vmin,
                np.max(campo) if vmax is None else vmax)
    ax.add_collection(pc)
    if malha:
        seg = [np.column_stack((ElemX[e, [0, 1, 2, 0]], ElemY[e, [0, 1, 2, 0]]))
               for e in range(nelem)]
        ax.add_collection(LineCollection(seg, colors="k", linewidths=0.2))
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    cb = plt.colorbar(pc, ax=ax, fraction=0.03, pad=0.02)
    if rotulo:
        cb.set_label(rotulo)
    if titulo:
        ax.set_title(titulo, fontsize=10)
    return pc


def plot_campo_nodal(SysCoor, ElemNode0, valores, ax=None, niveis=14,
                     cmap="viridis", titulo=None, rotulo=None, deslocado=None):
    """
    Contornos preenchidos de um campo nodal. Cada LST e dividido em 4
    triangulos lineares (usando os nos do meio dos lados).

    deslocado : (nnode, 2) opcional; coordenadas a usar no lugar de SysCoor
                (malha deformada).
    """
    ax = plt.gca() if ax is None else ax
    xy = SysCoor if deslocado is None else deslocado
    n = ElemNode0
    sub = np.vstack((n[:, [0, 3, 5]], n[:, [3, 1, 4]], n[:, [5, 4, 2]], n[:, [3, 4, 5]]))
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1], sub)
    cs = ax.tricontourf(tri, valores, levels=niveis, cmap=cmap)
    ax.set_aspect("equal", adjustable="box")
    cb = plt.colorbar(cs, ax=ax, fraction=0.03, pad=0.02)
    if rotulo:
        cb.set_label(rotulo)
    if titulo:
        ax.set_title(titulo, fontsize=10)
    return cs


def plot_direcoes_principais(ElemX, ElemY, sig, ax=None, passo=1, escala=None,
                             titulo=None):
    """
    Setas das direcoes principais de tensao no centroide de cada elemento
    (media dos pontos de Gauss): vermelho = s1 (menos compressiva),
    azul = s3 (mais compressiva). Comprimento proporcional a |tensao|.
    """
    ax = plt.gca() if ax is None else ax
    sm = sig.mean(axis=1)                               # (nelem, nsigma)
    inv = invariantes(sm)
    cx, cy = ElemX[:, :3].mean(axis=1), ElemY[:, :3].mean(axis=1)
    th = np.radians(inv["theta"])
    sel = np.arange(0, cx.size, passo)
    if escala is None:
        escala = 1.4 * np.sqrt(np.ptp(cx) * np.ptp(cy) / max(sel.size, 1)) \
            / max(np.max(np.abs(inv["s3"])), 1e-9)
    for s, cor, ang in ((inv["s1"], "tab:red", th), (inv["s3"], "tab:blue", th + np.pi / 2)):
        L = escala * np.abs(s[sel])
        dx, dy = L * np.cos(ang[sel]), L * np.sin(ang[sel])
        segs = [[(x - a, y - b), (x + a, y + b)] for x, y, a, b in
                zip(cx[sel], cy[sel], 0.5 * dx, 0.5 * dy)]
        ax.add_collection(LineCollection(segs, colors=cor, linewidths=1.2))
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    if titulo:
        ax.set_title(titulo, fontsize=10)
