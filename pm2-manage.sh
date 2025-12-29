#!/bin/bash

# PM2 Management Script for ssKIND Backend
# Usage: ./pm2-manage.sh [start|stop|restart|status|logs|delete]

APP_NAME="sskind-backend"
CONFIG_FILE="ecosystem.config.js"

case "$1" in
    start)
        echo "🚀 Starting ssKIND Backend with PM2..."
        pm2 start $CONFIG_FILE --env production
        pm2 save
        ;;
    start-dev)
        echo "🔧 Starting ssKIND Backend in development mode..."
        pm2 start $CONFIG_FILE --env development
        pm2 save
        ;;
    stop)
        echo "⏹️  Stopping ssKIND Backend..."
        pm2 stop $APP_NAME
        ;;
    restart)
        echo "🔄 Restarting ssKIND Backend..."
        pm2 restart $APP_NAME
        ;;
    reload)
        echo "🔄 Reloading ssKIND Backend (zero-downtime)..."
        pm2 reload $APP_NAME
        ;;
    status)
        echo "📊 ssKIND Backend Status:"
        pm2 status $APP_NAME
        ;;
    logs)
        echo "📋 Showing ssKIND Backend logs (Ctrl+C to exit)..."
        pm2 logs $APP_NAME
        ;;
    delete)
        echo "🗑️  Deleting ssKIND Backend from PM2..."
        pm2 delete $APP_NAME
        ;;
    list)
        echo "📋 All PM2 processes:"
        pm2 list
        ;;
    monit)
        echo "📊 Opening PM2 monitoring..."
        pm2 monit
        ;;
    *)
        echo "Usage: $0 {start|start-dev|stop|restart|reload|status|logs|delete|list|monit}"
        echo ""
        echo "Commands:"
        echo "  start      - Start in production mode (port 9117)"
        echo "  start-dev  - Start in development mode with reload"
        echo "  stop       - Stop the application"
        echo "  restart    - Restart the application"
        echo "  reload     - Zero-downtime reload"
        echo "  status     - Show application status"
        echo "  logs       - Show application logs"
        echo "  delete     - Remove from PM2"
        echo "  list       - List all PM2 processes"
        echo "  monit      - Open PM2 monitoring dashboard"
        exit 1
        ;;
esac
