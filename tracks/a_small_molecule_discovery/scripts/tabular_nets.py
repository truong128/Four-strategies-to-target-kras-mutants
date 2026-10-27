"""Four tabular neural nets for the docking-label benchmark.

Two FT-Transformer variants group fingerprint bits into tokens so attention
stays tractable at 2048 bits. RealMLP is a scaled deep MLP (SELU). Tabular
ResNet follows the residual-block design used for tabular data.
"""
import math
import numpy as np
import torch
from torch import nn
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.preprocessing import StandardScaler


def _device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


class _FTTransformer(nn.Module):
    def __init__(self, n_features, group, d_token, n_layers, n_heads, dropout):
        super().__init__()
        self.group = int(group)
        self.n_tokens = int(math.ceil(n_features / self.group))
        self.pad = self.n_tokens * self.group - int(n_features)
        self.weight = nn.Parameter(torch.empty(self.n_tokens, self.group, d_token))
        nn.init.xavier_uniform_(self.weight.view(self.n_tokens, -1))
        self.bias = nn.Parameter(torch.zeros(self.n_tokens, d_token))
        self.cls = nn.Parameter(torch.zeros(1, 1, d_token))
        layer = nn.TransformerEncoderLayer(
            d_model=d_token, nhead=n_heads, dim_feedforward=d_token * 4,
            dropout=dropout, batch_first=True, activation="gelu", norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, num_layers=n_layers)
        self.head = nn.Linear(d_token, 1)

    def forward(self, x):
        if self.pad:
            x = nn.functional.pad(x, (0, self.pad))
        b = x.shape[0]
        tok = torch.einsum("btg,tgd->btd", x.view(b, self.n_tokens, self.group), self.weight) + self.bias
        h = self.encoder(torch.cat([self.cls.expand(b, -1, -1), tok], dim=1))
        return self.head(h[:, 0]).squeeze(-1)


class _ResBlock(nn.Module):
    def __init__(self, width, dropout):
        super().__init__()
        self.norm = nn.BatchNorm1d(width)
        self.lin1 = nn.Linear(width, width)
        self.lin2 = nn.Linear(width, width)
        self.drop = nn.Dropout(dropout)

    def forward(self, x):
        h = self.drop(nn.functional.relu(self.lin1(self.norm(x))))
        h = self.drop(self.lin2(nn.functional.relu(h)))
        return x + h


class _TabularResNet(nn.Module):
    def __init__(self, n_features, width=128, n_blocks=3, dropout=0.1):
        super().__init__()
        self.proj = nn.Linear(n_features, width)
        self.blocks = nn.Sequential(*[_ResBlock(width, dropout) for _ in range(n_blocks)])
        self.head = nn.Linear(width, 1)

    def forward(self, x):
        return self.head(self.blocks(nn.functional.relu(self.proj(x)))).squeeze(-1)


class _RealMLPNet(nn.Module):
    def __init__(self, n_features, hidden=256, n_layers=3, dropout=0.1):
        super().__init__()
        layers = []
        d = n_features
        for _ in range(n_layers):
            lin = nn.Linear(d, hidden)
            nn.init.kaiming_normal_(lin.weight, nonlinearity="linear")
            layers += [lin, nn.SELU(), nn.AlphaDropout(dropout)]
            d = hidden
        self.backbone = nn.Sequential(*layers)
        self.head = nn.Linear(hidden, 1)

    def forward(self, x):
        return self.head(self.backbone(x)).squeeze(-1)


ARCH = {
    "ft_transformer": dict(group=32, d_token=64, n_layers=3, n_heads=8, dropout=0.1),
    "ft_transformer_lite": dict(group=64, d_token=32, n_layers=2, n_heads=4, dropout=0.1),
    "realmlp": dict(hidden=256, n_layers=3, dropout=0.1),
    "tabular_resnet": dict(width=128, n_blocks=3, dropout=0.1),
}


class TabularNetClassifier(BaseEstimator, ClassifierMixin):
    """sklearn-style wrapper. The network is moved to CPU after fit so joblib workers can score it."""

    def __init__(self, arch="realmlp", epochs=8, batch_size=2048, patience=3,
                 learning_rate=1e-3, max_rows=60000, seed=42):
        self.arch = arch
        self.epochs = epochs
        self.batch_size = batch_size
        self.patience = patience
        self.learning_rate = learning_rate
        self.max_rows = max_rows
        self.seed = seed

    def fit(self, X, y):
        if self.arch not in ARCH:
            raise ValueError(self.arch)
        X = np.asarray(X, dtype=np.float32)
        y = np.asarray(y).astype(np.int64).ravel()
        rng = np.random.default_rng(self.seed)
        if len(y) > int(self.max_rows):
            pos = np.where(y == 1)[0]
            neg = np.where(y == 0)[0]
            n_pos = min(len(pos), max(1, int(self.max_rows) // 2))
            take = np.concatenate([
                rng.choice(pos, n_pos, replace=False),
                rng.choice(neg, min(len(neg), int(self.max_rows) - n_pos), replace=False),
            ])
            rng.shuffle(take)
            X, y = X[take], y[take]
        self.classes_ = np.array([0, 1])
        self.scaler_ = StandardScaler().fit(X)
        Xs = self.scaler_.transform(X).astype(np.float32)
        n = len(y)
        perm = rng.permutation(n)
        n_val = max(1, int(0.15 * n)) if n >= 40 and len(np.unique(y)) == 2 else 0
        val_ix, tr_ix = perm[:n_val], perm[n_val:]
        if len(np.unique(y[tr_ix])) < 2:
            tr_ix, val_ix = perm, perm[:0]
        device = _device()
        torch.manual_seed(int(self.seed))
        spec = ARCH[self.arch]
        n_features = Xs.shape[1]
        if self.arch.startswith("ft_transformer"):
            net = _FTTransformer(n_features, **spec)
        elif self.arch == "tabular_resnet":
            net = _TabularResNet(n_features, **spec)
        else:
            net = _RealMLPNet(n_features, **spec)
        net = net.to(device)
        opt = torch.optim.AdamW(net.parameters(), lr=float(self.learning_rate), weight_decay=1e-4)
        n_pos = max(1, int((y[tr_ix] == 1).sum()))
        n_neg = max(1, int((y[tr_ix] == 0).sum()))
        loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([n_neg / n_pos], device=device))
        xt = torch.from_numpy(Xs[tr_ix])
        yt = torch.from_numpy(y[tr_ix].astype(np.float32))
        xv = torch.from_numpy(Xs[val_ix]) if len(val_ix) else None
        yv = y[val_ix] if len(val_ix) else None
        best_state, best_loss, wait = None, math.inf, 0
        bs = max(32, int(self.batch_size))
        net.train()
        for _ in range(int(self.epochs)):
            order = rng.permutation(len(tr_ix))
            for start in range(0, len(order), bs):
                b = order[start:start + bs]
                xb = xt[b].to(device, non_blocking=True)
                yb = yt[b].to(device, non_blocking=True)
                opt.zero_grad(set_to_none=True)
                loss = loss_fn(net(xb), yb)
                loss.backward()
                opt.step()
            if xv is None:
                continue
            net.eval()
            with torch.no_grad():
                logits = []
                for start in range(0, len(val_ix), bs):
                    logits.append(net(xv[start:start + bs].to(device)).float().cpu())
                prob = torch.sigmoid(torch.cat(logits)).numpy()
            # early stop on validation log loss
            p = np.clip(prob, 1e-6, 1 - 1e-6)
            vloss = float(-(yv * np.log(p) + (1 - yv) * np.log(1 - p)).mean())
            net.train()
            if vloss + 1e-4 < best_loss:
                best_loss, wait = vloss, 0
                best_state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}
            else:
                wait += 1
                if wait >= int(self.patience):
                    break
        if best_state is not None:
            net.load_state_dict(best_state)
        self.model_ = net.cpu().eval()
        return self

    def predict_proba(self, X):
        X = np.asarray(X, dtype=np.float32)
        xs = self.scaler_.transform(X).astype(np.float32)
        self.model_.eval()
        out = np.empty(len(xs), dtype=np.float32)
        bs = max(256, int(self.batch_size))
        with torch.no_grad():
            for start in range(0, len(xs), bs):
                logits = self.model_(torch.from_numpy(xs[start:start + bs]))
                out[start:start + bs] = torch.sigmoid(logits).numpy()
        out = np.clip(out, 0, 1)
        return np.column_stack([1 - out, out])

    def predict(self, X):
        return (self.predict_proba(X)[:, 1] >= 0.5).astype(int)
