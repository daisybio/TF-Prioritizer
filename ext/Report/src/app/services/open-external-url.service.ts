import {Injectable} from '@angular/core';

@Injectable({
  providedIn: 'root'
})
export class OpenExternalUrlService {

  constructor() {
  }

  openUrl(url: string) {
    window.open(url);
  }
}
