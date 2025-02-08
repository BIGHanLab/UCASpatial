import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import NMF
from sklearn.preprocessing import normalize
from scipy.optimize import nnls
from tqdm import tqdm
from scipy import stats
import os
import logging
import warnings
from typing import Optional, Union

logging.basicConfig(level=logging.INFO)
warnings.filterwarnings('ignore')

class UCASpatial
    def __init__(self, 
                 sc_ref ad.AnnData,
                 st_vis ad.AnnData,
                 clust_vr str,
                 assay str = 'RNA',
                 slot str = 'counts',
                 output_path Optional[str] = None,
                 cluster_markers Optional[pd.DataFrame] = None,
                 min_pct float = 0.2,
                 logfc_threshold float = 0.25,
                 min_diff_pct float = 0.1,
                 normalize_method str = 'uv',
                 downsample_n int = 1,
                 n_cluster int = 100,
                 n_top Optional[int] = None,
                 min_cont float = 0.01,
                 remove_RPL bool = False,
                 remove_MT bool = False,
                 cos_filter bool = True,
                 cos_mu float = 1,
                 cos_n_genes_user int = 900,
                 marker_slot str = 'data',
                 unit str = log2,
                 random_seed int = 10000,
                 meta_filter bool = True,
                 nmf_tol float = 1e-4,
                 meta_assay str = 'integrated',
                 meta_ndims int = 30,
                 meta_resolution int = 100,
                 meta_purity float = 0.95,
                 ent_filter_threshold float = 0.5,
                 cos_filter_threshold float = 0.05,
                 weight_filter_threshold float = 0.2)
        
        self.sc_ref = sc_ref
        self.st_vis = st_vis
        self.clust_vr = clust_vr
        self.assay = assay
        self.slot = slot
        self.output_path = output_path or os.getcwd()
        self.cluster_markers = cluster_markers
        self.min_pct = min_pct
        self.logfc_threshold = logfc_threshold
        self.min_diff_pct = min_diff_pct
        self.normalize_method = normalize_method
        self.downsample_n = downsample_n
        self.n_cluster = n_cluster
        self.n_top = n_top
        self.min_cont = min_cont
        self.remove_RPL = remove_RPL
        self.remove_MT = remove_MT
        self.cos_filter = cos_filter
        self.cos_mu = cos_mu
        self.cos_n_genes_user = cos_n_genes_user
        self.marker_slot = marker_slot
        self.unit = unit
        self.random_seed = random_seed
        self.meta_filter = meta_filter
        self.nmf_tol = nmf_tol
        self.meta_assay = meta_assay
        self.meta_ndims = meta_ndims
        self.meta_resolution = meta_resolution
        self.meta_purity = meta_purity
        self.ent_filter_threshold = ent_filter_threshold
        self.cos_filter_threshold = cos_filter_threshold
        self.weight_filter_threshold = weight_filter_threshold

    def preprocess(self)
        logging.info(Step0 Preprocessing...)
        self._validate_inputs()
        self._filter_genes()
        
    def _validate_inputs(self)
        if self.clust_vr not in self.sc_ref.obs
            raise ValueError(f{self.clust_vr} not found in sc_ref.obs)
            
        if not os.path.exists(self.output_path)
            os.makedirs(self.output_path)

    def _filter_genes(self)
        if self.remove_RPL
            rpl_genes = self.st_vis.var_names[self.st_vis.var_names.str.startswith(('RPL', 'RPS'))]
            self.st_vis = self.st_vis[, ~self.st_vis.var_names.isin(rpl_genes)].copy()
            
        if self.remove_MT
            mt_genes = self.st_vis.var_names[self.st_vis.var_names.str.startswith('MT-')]
            self.st_vis = self.st_vis[, ~self.st_vis.var_names.isin(mt_genes)].copy()

    def downsample_sc_ref(self)
        logging.info(Downsampling sc_ref...)
        np.random.seed(self.random_seed)
        
        clusters = self.sc_ref.obs[self.clust_vr].unique()
        sampled_cells = []
        
        for cluster in clusters
            cluster_cells = self.sc_ref.obs[self.sc_ref.obs[self.clust_vr] == cluster].index
            if len(cluster_cells) = self.n_cluster
                sampled_cells.extend(cluster_cells)
                continue
                
            min_variance = np.inf
            best_sample = None
            
            for _ in range(self.downsample_n)
                sample = np.random.choice(cluster_cells, size=self.n_cluster, replace=False)
                expr = self.sc_ref[sample, ].X
                
                if sp.issparse(expr)
                    expr = expr.toarray()
            
                variance = expr.var(axis=0).sum()
                
                if variance  min_variance
                    min_variance = variance
                    best_sample = sample
            
            sampled_cells.extend(best_sample)
        
        self.sc_ref = self.sc_ref[sampled_cells].copy()

    def find_markers(self)
        logging.info(Finding marker genes...)
        
        self.sc_ref.raw = self.sc_ref
        
        self.sc_ref.X = self.sc_ref.layers[self.marker_slot]
        
        sc.pp.normalize_total(self.sc_ref, target_sum=1e4)
        sc.pp.log1p(self.sc_ref)
        
        sc.tl.rank_genes_groups(
            self.sc_ref,
            groupby=self.clust_vr,
            method='wilcoxon',
            n_genes=self.sc_ref.shape[1],
            use_raw=False,
            layer=self.marker_slot,
            key_added='dea',
            pts=True
        )
        
        markers = []
        for cluster in self.sc_ref.obs[self.clust_vr].unique()
            df = sc.get.rank_genes_groups_df(self.sc_ref, 
                                       group=cluster,
                                       key='dea',
                                       pval_cutoff=0.05,
                                       log2fc_min=self.logfc_threshold)
            df = df.rename(columns={
                'names' 'gene',
                'logfoldchanges' 'avg_logFC',
                'pct_nz_group' 'pct.1',
                'pct_nz_reference' 'pct.2',
                'scores' 'score',
                'pvals_adj' 'p_val_adj'
            })
            df['cluster'] = cluster
            markers.append(df)
            
        self.cluster_markers = pd.concat(markers)
        self.cluster_markers = self.cluster_markers.sort_values(by=['gene', 'score'], ascending=[True, False]).drop_duplicates('gene')
        self.cluster_markers.reset_index(drop=True, inplace=True)
        self._filter_markers() 
#         print(self.sc_ref.layers.keys())
        self._calculate_weights()

    def _filter_markers(self)
        self.cluster_markers = self.cluster_markers[
            (self.cluster_markers['p_val_adj']  0.05) &
            (self.cluster_markers['avg_logFC']  self.logfc_threshold) &
            (self.cluster_markers['pct.1']  self.min_pct) &
            (self.cluster_markers['pct.1'] - self.cluster_markers['pct.2']  self.min_diff_pct)
        ]
    
    def _calculate_weights(self)
        logging.info(Calculating entropy weights...)
        # Simplified entropy calculation
        clusters = self.sc_ref.obs[self.clust_vr].unique()
        
        # 仅保留cluster_markers中的基因
        var_names = self.sc_ref.var_names
        genes_of_interest = self.cluster_markers['gene'].unique()
        indices_of_interest = var_names.get_indexer(genes_of_interest)
        var_names = var_names[indices_of_interest]
        
        expr = self.sc_ref.layers[self.marker_slot]
        
        if sp.issparse(expr)
            expr = expr.toarray()
        expr = expr[, indices_of_interest]
        
        
        pct = pd.DataFrame(index=var_names)
        for cluster in clusters
            cluster_cells = self.sc_ref.obs[self.sc_ref.obs[self.clust_vr] == cluster].index
            cluster_cells_idx = self.sc_ref.obs.index.get_indexer(cluster_cells)
            pct[cluster] = (expr[cluster_cells_idx,]  0).mean(axis=0)
        entropy = stats.entropy(pct.T, base=2)
        max_entropy = np.log2(len(clusters))
        entropy_adj = max_entropy - entropy
        
        self.cluster_markers['ent_adj'] = entropy_adj
        self.cluster_markers['weight'] = (
            self.cluster_markers['ent_adj'] 
            self.cluster_markers['pct.1'] 
            self.cluster_markers['avg_logFC']
        )

    def train_nmf(self)
        logging.info(Training NMF model...)
        common_genes = list(set(self.sc_ref.var_names) & set(self.st_vis.var_names) & set(self.cluster_markers['gene']))
        self.sc_ref = self.sc_ref[, common_genes].copy()
        self.st_vis = self.st_vis[, common_genes].copy()
        
        if self.normalize_method == 'uv'
            expr = normalize(self.sc_ref.X, norm='l2', axis=0)
        else
            expr = self.sc_ref.X
            
        self.nmf = NMF(
            n_components=len(self.sc_ref.obs[self.clust_vr].unique()),
            tol=self.nmf_tol,
            random_state=self.random_seed
        )
        self.W = self.nmf.fit_transform(expr)
        self.H = self.nmf.components_

#     def deconvolute(self)
#         logging.info(Performing deconvolution...)
#         if self.normalize_method == 'uv'
#             st_expr = normalize(self.st_vis.X, norm='l2', axis=0)
#         else
#             st_expr = self.st_vis.X
            
#         C = cosine_similarity(self.H.T, st_expr.T)
#         proportions = np.zeros((st_expr.shape[0], self.H.shape[0]))
        
#         for i in range(st_expr.shape[0])
#             proportions[i] = nnls(self.W, st_expr[i])[0]
            
#         proportions = proportions  proportions.sum(axis=1, keepdims=True)
#         proportions[proportions  self.min_cont] = 0
#         proportions = proportions  proportions.sum(axis=1, keepdims=True)
        
#         self.proportions = pd.DataFrame(
#             proportions,
#             index=self.st_vis.obs_names,
#             columns=self.sc_ref.obs[self.clust_vr].unique()
#         )
        
    def deconvolute(self)
        logging.info(Performing weighted NNLS deconvolution...)

        # 获取公共基因
        common_genes = list(set(self.sc_ref.var_names) & set(self.st_vis.var_names) & set(self.cluster_markers['gene']))
        st_expr = self.st_vis[, common_genes].X.T  # 转换为 (genes, spots)

        # Step 1 标准化空间表达矩阵
        if self.normalize_method == 'cpm'
            count_matr = self._cpm_normalize(st_expr)
        elif self.normalize_method == 'uv'
            count_matr = self._unit_variance_normalize(st_expr)
        elif self.normalize_method == 'raw'
            count_matr = st_expr
        else
            raise ValueError(fUnsupported normalize method {self.normalize_method})

        # Step 2 构建权重矩阵
        weight_df = (
            self.cluster_markers
            .groupby('gene')
            .apply(lambda x x.nlargest(1, 'weight'))
            .reset_index(drop=True)[['gene', 'weight']]
            .drop_duplicates()
        )
        weight_df = weight_df[weight_df['gene'].isin(common_genes)]

        # 处理无穷大权重
        max_weight = weight_df.loc[weight_df['weight'] != np.inf, 'weight'].max()
        weight_df.loc[weight_df['weight'] == np.inf, 'weight'] = max_weight

        # 创建对角权重矩阵
        weight_vec = np.sqrt(weight_df.set_index('gene').reindex(common_genes)['weight'].values)
        weight_mat = np.diag(weight_vec)

        # Step 3 调整NMF矩阵和表达矩阵
        W = self.nmf.components_.T  # (genes, topics)
        W_weighted = weight_mat @ W
        count_weighted = weight_mat @ count_matr.T

        # Step 4 NNLS求解
        n_topics = W_weighted.shape[1]
        n_spots = count_weighted.shape[1]

        coef_pred = np.zeros((n_topics, n_spots))
        for i in tqdm(range(n_spots), desc=Deconvoluting spots)
            coef, _ = nnls(W_weighted, count_weighted[, i])
            coef_pred[, i] = coef

        # Step 5 计算最终比例
        reference_profiles = self._get_topic_profiles()
        decon_res = self._calculate_proportions(coef_pred, reference_profiles)

        self.proportions = decon_res
        return decon_res

    def _cpm_normalize(self, expr)
        CPM标准化
        lib_size = expr.sum(axis=0)
        return expr  lib_size  1e6

    def _unit_variance_normalize(self, expr)
        单位方差标准化
        if isinstance(expr, np.ndarray)
            stds = np.std(expr, axis=0)
        elif isinstance(expr, sp.spmatrix)
            stds = np.sqrt(np.mean(expr.power(2), axis=0))
        else
            raise ValueError(Unsupported type for expr)
#         stds = np.std(expr, axis=1, keepdims=True)
        stds[stds == 0] = 1  # 避免除以0
        return (expr  stds).T  # 转置回 (genes, spots)

    def _get_topic_profiles(self)
        获取主题-细胞类型特征矩阵
        # 假设已通过NMF得到H矩阵 (topics x cells)
        H = self.nmf.transform(self.sc_ref.X)  # 需要根据实际NMF实现调整

        # 按细胞类型聚合
        cell_types = self.sc_ref.obs[self.clust_vr].values
        unique_types = np.unique(cell_types)
        profiles = np.zeros((len(unique_types), H.shape[1]))

        for i, ct in enumerate(unique_types)
            profiles[i] = H[cell_types == ct].mean(axis=0)

        return profiles  # (cell_types, topics)

    def _calculate_proportions(self, coef_pred, reference_profiles)
        计算最终细胞类型比例
        n_spots = coef_pred.shape[1]
        n_celltypes = reference_profiles.shape[0]

        proportions = np.zeros((n_spots, n_celltypes + 1))  # 最后列为残差

        for i in tqdm(range(n_spots), desc=Calculating proportions)
            # 计算细胞类型组成
            weights, _ = nnls(reference_profiles.T, coef_pred[, i])
            comp = weights  weights.sum()

            # 应用最小贡献阈值
            comp[comp  self.min_cont] = 0
            comp = comp  comp.sum()

            # 计算残差
            pred = reference_profiles.T @ comp
            total_ss = np.sum((coef_pred[, i])2)
            res_ss = np.sum((coef_pred[, i] - pred)2)  total_ss

            proportions[i, -1] = comp
            proportions[i, -1] = res_ss

        # 转换为DataFrame
        cell_types = np.unique(self.sc_ref.obs[self.clust_vr])
        return pd.DataFrame(
            proportions,
            index=self.st_vis.obs_names,
            columns=list(cell_types) + ['residual_ss']
        )

    def run(self)
        self.preprocess()
        
        if self.meta_filter
            self._meta_filter()
            
        if self.downsample_n  0
            self.downsample_sc_ref()
            
        if self.cluster_markers is None
            self.find_markers()
            
        self.train_nmf()
        self.deconvolute()
        
        return {
            'nmf_components' self.H,
            'proportions' self.proportions
        }

    def _meta_filter(self)
        logging.info(Performing meta filtering...)
        sc_ref = self.sc_ref.copy()
        
        # Step 1 设置当前assay
        current_assay = self.meta_assay
        if current_assay not in sc_ref.layers
            sc_ref.X = sc_ref.layers[current_assay] = sc_ref.X.copy()
        
        # Step 2 检查降维数据是否存在
        if 'X_pca' not in sc_ref.obsm or f'{current_assay}_pca' not in sc_ref.uns
            logging.warning(No PCA found, performing PCA...)
            sc.pp.highly_variable_genes(
                sc_ref, 
                flavor='seurat_v3',
                n_top_genes=2000,
                layer=current_assay
            )
            sc.pp.scale(sc_ref, layer=current_assay)
            sc.tl.pca(
                sc_ref, 
                n_comps=50,
                use_highly_variable=True,
                layer=current_assay,
                svd_solver='arpack'
            )
        # Step 3 检查邻居图是否存在
        if 'neighbors' not in sc_ref.uns
            logging.warning(No neighbor graph found, constructing...)
            sc.pp.neighbors(
                sc_ref,
                n_pcs=self.meta_ndims,
                use_rep='X_pca',
                random_state=self.random_seed
            )
        # Step 4 进行Leiden聚类
        cluster_key = f{current_assay}_leiden_{self.meta_resolution}
        if cluster_key not in sc_ref.obs
            sc.tl.leiden(
                sc_ref,
                resolution=self.meta_resolution,
                key_added=cluster_key,
                random_state=self.random_seed
            )
        # 创建Meta_cell列
        meta_cell_col = fMeta_cell_{self.meta_resolution}
        sc_ref.obs[meta_cell_col] = [
            fMeta_cell_{x} 
            for x in sc_ref.obs[cluster_key]
        ]
        # Step 5 计算元细胞纯度
        cross_tab = pd.crosstab(
            sc_ref.obs[self.clust_vr],
            sc_ref.obs[meta_cell_col]
        )
        
        # 计算每个Meta_cell的纯度
        purity_df = cross_tab.apply(lambda x x.max()  x.sum(), axis=0).reset_index()
        purity_df.columns = ['Meta_cell', 'Purity']
        
        # 过滤纯度达标的Meta_cell
        valid_meta_cells = purity_df[purity_df['Purity'] = self.meta_purity]['Meta_cell'].tolist()
        
        # Step 6 保留属于有效Meta_cell的细胞
        filtered_cells = sc_ref.obs[
            sc_ref.obs[meta_cell_col].isin(valid_meta_cells)
        ].index
    
        # 更新sc_ref并恢复原始assay
        self.sc_ref = sc_ref[filtered_cells].copy()
        self.sc_ref.X = self.sc_ref.layers[self.meta_assay]  # 恢复原始数据

        # 保存结果（可选）
        if self.output_path
            import os
            output_file = os.path.join(self.output_path, 'sc_ref_meta_filter.h5ad')
            self.sc_ref.write(output_file)
            logging.info(fFiltered data saved to {output_file})