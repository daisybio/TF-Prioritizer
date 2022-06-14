import { ComponentFixture, TestBed } from '@angular/core/testing';

import { CrossLinksComponent } from './cross-links.component';

describe('CrossLinksComponent', () => {
  let component: CrossLinksComponent;
  let fixture: ComponentFixture<CrossLinksComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ CrossLinksComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(CrossLinksComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
