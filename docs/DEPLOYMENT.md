# HypoMap Backend Deployment Guide

This guide covers deploying the HypoMap backend on a production server using micromamba and PM2.

## Prerequisites

- Linux server (Ubuntu 20.04+ recommended)
- Node.js 18+ (for PM2)
- Git
- Nginx (for reverse proxy)

## 1. Install Micromamba

```bash
# Download and install micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Move to a permanent location
sudo mv bin/micromamba /usr/local/bin/

# Initialize micromamba
micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc

# Verify installation
micromamba --version
```

## 2. Create Conda Environment

```bash
# Create environment with Python 3.11
micromamba create -n hypomap python=3.11 -c conda-forge -y

# Activate environment
micromamba activate hypomap

# Verify Python
python --version
which python
# Should show: ~/micromamba/envs/hypomap/bin/python
```

## 3. Clone and Setup Repository

```bash
# Clone repository
cd /opt  # or your preferred location
sudo git clone https://github.com/your-org/HypoMap.git
sudo chown -R $USER:$USER HypoMap

# Navigate to backend
cd HypoMap/hypomap-backend

# Install Python dependencies
micromamba activate hypomap
pip install -r requirements.txt
```

## 4. Configure Environment

```bash
# Copy and edit environment file
cp env.example .env
nano .env
```

Edit `.env` for production:

```env
# HypoMap Backend Configuration
PORT=9120
HOST=0.0.0.0
DEBUG=False
RELOAD=False

# Database Configuration
DATABASE_URL=sqlite:///./app.db

# API Configuration
GLOBAL_PREFIX=/hypomap-backend
API_PREFIX=/api/v1
DOCS_URL=/docs
REDOC_URL=/redoc

# CORS Configuration
CORS_ORIGINS=["https://your-domain.com"]
CORS_CREDENTIALS=True
CORS_METHODS=["*"]
CORS_HEADERS=["*"]

# App Metadata
APP_TITLE=HypoMap Backend API
APP_DESCRIPTION=Single-cell RNA-seq analysis with AI-powered hypothesis generation
APP_VERSION=1.0.0
```

## 5. Setup PM2 for Process Management

```bash
# Install PM2 globally
npm install -g pm2

# Create PM2 ecosystem config
nano ecosystem.config.js
```

Create `ecosystem.config.js`:

```javascript
module.exports = {
  apps: [
    {
      name: "hypomap-backend",
      script: "main.py",
      interpreter: "/home/YOUR_USER/micromamba/envs/hypomap/bin/python",
      cwd: "/opt/HypoMap/hypomap-backend",
      env: {
        HOST: "0.0.0.0",
        PORT: "9120",
        NODE_ENV: "production",
        DEBUG: "False",
        RELOAD: "False",
        DATABASE_URL: "sqlite:///./app.db",
        GLOBAL_PREFIX: "/hypomap-backend",
        API_PREFIX: "/api/v1",
        DOCS_URL: "/docs",
        REDOC_URL: "/redoc",
        CORS_ORIGINS: '["*"]',
        CORS_CREDENTIALS: "True",
        CORS_METHODS: '["*"]',
        CORS_HEADERS: '["*"]',
        APP_TITLE: "HypoMap Backend API",
        APP_DESCRIPTION: "Single-cell RNA-seq analysis with AI",
        APP_VERSION: "1.0.0",
      },
      instances: 1,
      exec_mode: "fork",
      watch: false,
      max_memory_restart: "2G",
      error_file: "./logs/err.log",
      out_file: "./logs/out.log",
      log_file: "./logs/combined.log",
      time: true,
      autorestart: true,
      restart_delay: 1000,
      max_restarts: 10,
      min_uptime: "10s",
      kill_timeout: 5000,
      merge_logs: true,
      log_date_format: "YYYY-MM-DD HH:mm:ss Z",
    },
  ],
};
```

**Important:** Update the `interpreter` path to match your micromamba environment:
```bash
# Find your Python path
micromamba activate hypomap
which python
# Use this path in the interpreter field
```

## 6. Setup Log Directory

```bash
mkdir -p logs
chmod 755 logs
```

## 7. Prepare H5AD Data

```bash
# Create data directories
mkdir -p h5ad/raw h5ad/precomputed

# Copy your h5ad files
cp /path/to/your/data.h5ad h5ad/raw/

# Preprocess the data (generates UMAP, metadata, etc.)
micromamba activate hypomap
python scripts/preprocess_h5ad.py h5ad/raw/your_dataset.h5ad
```

## 8. Start the Service

```bash
# Start with PM2
pm2 start ecosystem.config.js

# Check status
pm2 status

# View logs
pm2 logs hypomap-backend

# Save PM2 process list (for auto-restart on reboot)
pm2 save

# Setup PM2 to start on boot
pm2 startup
# Follow the instructions printed by this command
```

## 9. Configure Nginx Reverse Proxy

### Option A: Adding to Existing BMBLX Server Configuration

If deploying to the bmblx.bmi.osumc.edu server, add the HypoMap locations to the existing nginx config:

```bash
sudo nano /etc/nginx/sites-available/bmblx
```

Add these location blocks inside the existing `server { ... }` block:

```nginx
    # ============================================
    # HypoMap - Single-cell Analysis with AI
    # ============================================

    # HypoMap Frontend
    location /hypomap {
        proxy_pass http://127.0.0.1:9121/hypomap;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # WebSocket support for Vite HMR (dev) or general WebSocket needs
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }

    # HypoMap Backend API
    location /hypomap-backend {
        proxy_pass http://127.0.0.1:9120/hypomap-backend;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # Increased timeouts for AI/LLM operations
        proxy_connect_timeout 600s;
        proxy_send_timeout 600s;
        proxy_read_timeout 600s;

        # For large file uploads (h5ad files)
        client_max_body_size 1G;
    }
```

### Option B: Full BMBLX-style Nginx Configuration

For reference, here's the complete SSL server block structure used on bmblx:

```nginx
# Connection upgrade map (add at top of nginx.conf or in http block)
map $http_upgrade $connection_upgrade {
    default upgrade;
    '' close;
}

server {
    listen 443 http2 ssl;
    listen [::]:443 http2 ssl;

    server_name bmblx.bmi.osumc.edu;

    # Self-signed SSL certificates
    ssl_certificate /etc/ssl/certs/nginx-selfsigned.crt;
    ssl_certificate_key /etc/ssl/private/nginx-selfsigned.key;
    ssl_dhparam /etc/ssl/certs/dhparam.pem;

    # Allow large file uploads
    client_max_body_size 1G;

    # ... other existing locations ...

    # ============================================
    # HypoMap - Single-cell Analysis with AI
    # ============================================

    # HypoMap Frontend (port 9121)
    location /hypomap {
        proxy_pass http://127.0.0.1:9121/hypomap;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # WebSocket support
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }

    # HypoMap Backend API (port 9120)
    location /hypomap-backend {
        proxy_pass http://127.0.0.1:9120/hypomap-backend;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # Increased timeouts for AI/LLM operations
        proxy_connect_timeout 600s;
        proxy_send_timeout 600s;
        proxy_read_timeout 600s;

        # Large file uploads
        client_max_body_size 1G;
    }

    # ... other existing locations ...
}
```

### Generate Self-Signed SSL Certificates (if needed)

```bash
# Generate SSL certificate
sudo openssl req -x509 -nodes -days 365 -newkey rsa:2048 \
    -keyout /etc/ssl/private/nginx-selfsigned.key \
    -out /etc/ssl/certs/nginx-selfsigned.crt

# Generate DH parameters (takes a few minutes)
sudo openssl dhparam -out /etc/ssl/certs/dhparam.pem 2048
```

### Test and Reload Nginx

```bash
# Test configuration syntax
sudo nginx -t

# If test passes, reload nginx
sudo systemctl reload nginx

# Or restart if needed
sudo systemctl restart nginx

# Check status
sudo systemctl status nginx
```

### Verify Deployment

```bash
# Test backend API
curl -k https://bmblx.bmi.osumc.edu/hypomap-backend/api/v1/health

# Test AI status
curl -k https://bmblx.bmi.osumc.edu/hypomap-backend/api/v1/ai/status

# Test frontend (should return HTML)
curl -k https://bmblx.bmi.osumc.edu/hypomap/ | head -20
```

## 10. PM2 Management Commands

```bash
# View all processes
pm2 list

# View logs
pm2 logs hypomap-backend
pm2 logs hypomap-backend --lines 100

# Restart
pm2 restart hypomap-backend

# Stop
pm2 stop hypomap-backend

# Delete from PM2
pm2 delete hypomap-backend

# Monitor resources
pm2 monit
```

## 11. Updating the Application

```bash
cd /opt/HypoMap/hypomap-backend

# Pull latest code
git pull origin main

# Activate environment and update dependencies
micromamba activate hypomap
pip install -r requirements.txt

# Restart the service
pm2 restart hypomap-backend

# Check logs for errors
pm2 logs hypomap-backend --lines 50
```

## 12. Ollama AI Service (Optional)

If using local Ollama instead of the bmblx server:

```bash
# Install Ollama
curl -fsSL https://ollama.com/install.sh | sh

# Pull required models
ollama pull qwen3:30b
ollama pull gpt-oss:20b
ollama pull qwen3:0.6b

# Start Ollama service
sudo systemctl enable ollama
sudo systemctl start ollama
```

Update `app/services/ollama_service.py` to point to local Ollama:
```python
OLLAMA_BASE_URL = "http://localhost:11434"
```

## Troubleshooting

### Check if backend is running
```bash
curl http://localhost:9120/hypomap-backend/api/v1/health
```

### Check PM2 logs for errors
```bash
pm2 logs hypomap-backend --err --lines 100
```

### Verify Python environment
```bash
micromamba activate hypomap
python -c "import fastapi; import anndata; print('OK')"
```

### Check port is not in use
```bash
lsof -i :9120
```

### Restart everything
```bash
pm2 restart all
sudo systemctl restart nginx
```

## Environment File Reference

Full `.env` template:

```env
# Server
PORT=9120
HOST=0.0.0.0
DEBUG=False
RELOAD=False

# Database
DATABASE_URL=sqlite:///./app.db

# API Routes
GLOBAL_PREFIX=/hypomap-backend
API_PREFIX=/api/v1
DOCS_URL=/docs
REDOC_URL=/redoc

# CORS (adjust for production)
CORS_ORIGINS=["https://your-domain.com"]
CORS_CREDENTIALS=True
CORS_METHODS=["*"]
CORS_HEADERS=["*"]

# App Info
APP_TITLE=HypoMap Backend API
APP_DESCRIPTION=Single-cell RNA-seq analysis with AI-powered hypothesis generation
APP_VERSION=1.0.0

# Ollama AI Service
OLLAMA_BASE_URL=https://bmblx.bmi.osumc.edu/ollama
```

## Security Checklist

- [ ] Set `DEBUG=False` in production
- [ ] Configure CORS_ORIGINS to only allow your frontend domain
- [ ] Use HTTPS with valid SSL certificates
- [ ] Set up firewall (ufw) to only allow necessary ports
- [ ] Regular security updates: `sudo apt update && sudo apt upgrade`
- [ ] Monitor logs for suspicious activity
